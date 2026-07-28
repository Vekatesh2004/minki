"""
Drug Matcher Module - Phase 4
Handles drug-target matching using DrugBank API and PharmGKB data
"""

import asyncio
import logging
from typing import Any, Dict, List, Optional, Tuple
import aiohttp
import json
import sqlite3
from pathlib import Path

logger = logging.getLogger(__name__)

class DrugMatcher:
    def __init__(self, config: Dict[str, Any]):
        """Initialize drug matcher with configuration"""
        self.config = config
        self.drugbank_config = config.get("drugbank", {})
        self.pharmgkb_config = config.get("pharmgkb", {})
        
        self.drugbank_api_url = self.drugbank_config.get("api_url", "https://go.drugbank.com/api/v1")
        self.drugbank_api_key = self.drugbank_config.get("api_key")
        self.rate_limit_delay = self.drugbank_config.get("rate_limit_delay", 1.0)
        
        self.pharmgkb_base_url = self.pharmgkb_config.get("base_url", "https://api.clinpgx.org/v1")
        self._pharmgkb_last_status = "not_checked"
        
        # Initialize the local public pharmacogene index. INSERT OR IGNORE plus
        # a uniqueness index keeps startup idempotent on fresh installations.
        self.db_path = str(Path(__file__).resolve().parents[1] / "drug_cache.db")
        self._init_local_database()
        
    def _init_local_database(self):
        """Initialize SQLite database for caching drug information"""
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Create tables for drug-target relationships
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS drug_targets (
                    id INTEGER PRIMARY KEY,
                    drug_id TEXT,
                    drug_name TEXT,
                    target_gene_symbol TEXT,
                    target_uniprot_id TEXT,
                    action TEXT,
                    organism TEXT,
                    pharmacology TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # Create table for pharmacogenomic annotations
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS pharmacogenomic_annotations (
                    id INTEGER PRIMARY KEY,
                    gene_symbol TEXT,
                    variant_hgvs TEXT,
                    drug_name TEXT,
                    phenotype TEXT,
                    evidence_level TEXT,
                    source TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # Remove legacy duplicate seed rows, then enforce idempotent seeds.
            cursor.execute('''
                DELETE FROM drug_targets
                WHERE id NOT IN (
                    SELECT MIN(id) FROM drug_targets
                    GROUP BY drug_id, drug_name, target_gene_symbol
                )
            ''')
            cursor.execute('''
                CREATE UNIQUE INDEX IF NOT EXISTS idx_unique_drug_target
                ON drug_targets(drug_id, drug_name, target_gene_symbol)
            ''')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_gene_symbol ON drug_targets(target_gene_symbol)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_uniprot_id ON drug_targets(target_uniprot_id)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_pgx_gene ON pharmacogenomic_annotations(gene_symbol)')
            self._seed_known_pharmacogenes(cursor)
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logger.error(f"Error initializing drug database: {str(e)}")

    @staticmethod
    def _seed_known_pharmacogenes(cursor):
        """Install the small public demonstration index on every fresh DB."""
        known = {
            "CYP2D6": [("Codeine", "substrate"), ("Tramadol", "substrate"),
                        ("Metoprolol", "substrate"), ("Risperidone", "substrate")],
            "CYP2C19": [("Clopidogrel", "substrate"), ("Omeprazole", "substrate"),
                         ("Escitalopram", "substrate")],
            "CYP2C9": [("Warfarin", "substrate"), ("Phenytoin", "substrate"),
                        ("Celecoxib", "substrate")],
            "TPMT": [("Azathioprine", "substrate"), ("Mercaptopurine", "substrate")],
            "DPYD": [("Fluorouracil", "substrate"), ("Capecitabine", "substrate")],
            "VKORC1": [("Warfarin", "target")],
        }
        rows = []
        for gene, drugs in known.items():
            for drug_name, action in drugs:
                rows.append((
                    f"pgx_{drug_name.lower().replace(' ', '_')}", drug_name, gene,
                    None, action, "Human", "Vitamin K recycling" if gene == "VKORC1" else "Metabolism",
                ))
        cursor.executemany('''
            INSERT OR IGNORE INTO drug_targets
                (drug_id, drug_name, target_gene_symbol, target_uniprot_id,
                 action, organism, pharmacology)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', rows)
            
    async def match_drugs(self, gene_symbol: str, uniprot_id: Optional[str] = None, 
                         hgvsp: Optional[str] = None) -> Dict[str, Any]:
        """
        Find drugs that target the specified gene/protein
        
        Args:
            gene_symbol: Gene symbol (e.g., CYP2D6)
            uniprot_id: Optional UniProt ID for additional matching
            hgvsp: Optional HGVSp notation for variant-specific lookups
            
        Returns:
            Dictionary containing matched drugs and pharmacogenomic information
        """
        try:
            # Search local cache first
            cached_drugs = await self._search_cached_drugs(gene_symbol, uniprot_id)
            
            # Search DrugBank API for additional matches
            drugbank_drugs = await self._search_drugbank_api(gene_symbol, uniprot_id)
            
            # Search ClinPGx/PharmGKB for pharmacogenomic annotations.
            # Availability is tracked separately from an empty successful result.
            pharmgkb_annotations = await self._search_pharmgkb(gene_symbol, hgvsp)
            pharmgkb_status = self._pharmgkb_last_status
            
            # Combine and rank results
            all_drugs = self._combine_drug_results(cached_drugs, drugbank_drugs)
            ranked_drugs = self._rank_drug_matches(
                all_drugs, pharmgkb_annotations, hgvsp,
                evidence_available=pharmgkb_status == "available",
            )
            
            return {
                "gene_symbol": gene_symbol,
                "uniprot_id": uniprot_id,
                "hgvsp": hgvsp,
                "matched_drugs": ranked_drugs,
                "pharmgkb_annotations": pharmgkb_annotations,
                "pharmgkb_status": pharmgkb_status,
                "total_matches": len(ranked_drugs)
            }
            
        except Exception as e:
            logger.error(f"Error matching drugs for gene {gene_symbol}: {str(e)}")
            return {
                "gene_symbol": gene_symbol,
                "uniprot_id": uniprot_id,
                "hgvsp": hgvsp,
                "error": str(e),
                "matched_drugs": [],
                "pharmgkb_annotations": [],
                "total_matches": 0
            }
            
    async def _search_cached_drugs(self, gene_symbol: str, uniprot_id: Optional[str]) -> List[Dict[str, Any]]:
        """Search for drugs in local cache"""
        
        cached_drugs = []
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            # Search by gene symbol
            cursor.execute('''
                SELECT drug_id, drug_name, target_gene_symbol, target_uniprot_id, 
                       action, organism, pharmacology
                FROM drug_targets 
                WHERE target_gene_symbol = ? OR target_gene_symbol LIKE ?
            ''', (gene_symbol, f'%{gene_symbol}%'))
            
            results = cursor.fetchall()
            
            for result in results:
                drug_info = {
                    "drug_id": result[0],
                    "drug_name": result[1],
                    "target_gene_symbol": result[2],
                    "target_uniprot_id": result[3],
                    "action": result[4],
                    "organism": result[5],
                    "pharmacology": result[6],
                    "source": "cache"
                }
                cached_drugs.append(drug_info)
                
            # Search by UniProt ID if provided
            if uniprot_id:
                cursor.execute('''
                    SELECT drug_id, drug_name, target_gene_symbol, target_uniprot_id, 
                           action, organism, pharmacology
                    FROM drug_targets 
                    WHERE target_uniprot_id = ?
                ''', (uniprot_id,))
                
                uniprot_results = cursor.fetchall()
                for result in uniprot_results:
                    drug_info = {
                        "drug_id": result[0],
                        "drug_name": result[1],
                        "target_gene_symbol": result[2],
                        "target_uniprot_id": result[3],
                        "action": result[4],
                        "organism": result[5],
                        "pharmacology": result[6],
                        "source": "cache"
                    }
                    cached_drugs.append(drug_info)
                    
            conn.close()
            
        except Exception as e:
            logger.error(f"Error searching cached drugs: {str(e)}")
            
        return cached_drugs
        
    async def _search_drugbank_api(self, gene_symbol: str, uniprot_id: Optional[str]) -> List[Dict[str, Any]]:
        """Search DrugBank API for drug-target relationships"""
        
        if not self.drugbank_api_key or self.drugbank_api_key == "your_api_key_here":
            logger.warning("DrugBank API key not configured, skipping API search")
            return []
            
        drugbank_drugs = []
        
        try:
            headers = {
                'Authorization': f'Bearer {self.drugbank_api_key}',
                'Accept': 'application/json'
            }
            
            # Search by gene symbol
            search_url = f"{self.drugbank_api_url}/targets"
            params = {'q': gene_symbol}
            
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=30)) as session:
                async with session.get(search_url, headers=headers, params=params) as response:
                    
                    if response.status == 200:
                        data = await response.json()
                        targets = data.get('targets', [])
                        
                        for target in targets:
                            # Get associated drugs for this target
                            target_id = target.get('id')
                            if target_id:
                                drugs = await self._get_drugs_for_target(session, headers, target_id)
                                drugbank_drugs.extend(drugs)
                                
                    elif response.status == 429:
                        logger.warning("DrugBank API rate limit hit")
                        await asyncio.sleep(self.rate_limit_delay * 2)
                    else:
                        logger.warning(f"DrugBank API error: {response.status}")
                        
        except Exception as e:
            logger.error(f"Error searching DrugBank API: {str(e)}")
            
        return drugbank_drugs
        
    async def _get_drugs_for_target(self, session: aiohttp.ClientSession, 
                                   headers: Dict[str, str], target_id: str) -> List[Dict[str, Any]]:
        """Get drugs associated with a specific target"""
        
        drugs = []
        
        try:
            drugs_url = f"{self.drugbank_api_url}/targets/{target_id}/drugs"
            
            async with session.get(drugs_url, headers=headers) as response:
                if response.status == 200:
                    data = await response.json()
                    drug_list = data.get('drugs', [])
                    
                    for drug in drug_list:
                        drug_info = {
                            "drug_id": drug.get('id'),
                            "drug_name": drug.get('name'),
                            "target_id": target_id,
                            "action": drug.get('action'),
                            "pharmacology": drug.get('general_function'),
                            "source": "drugbank_api"
                        }
                        drugs.append(drug_info)
                        
            await asyncio.sleep(self.rate_limit_delay)  # Rate limiting
            
        except Exception as e:
            logger.error(f"Error getting drugs for target {target_id}: {str(e)}")
            
        return drugs
        
    async def _search_pharmgkb(self, gene_symbol: str, hgvsp: Optional[str]) -> List[Dict[str, Any]]:
        """Search ClinPGx and distinguish an empty result from an outage."""

        annotations: List[Dict[str, Any]] = []
        self._pharmgkb_last_status = "not_checked"
        search_url = f"{self.pharmgkb_base_url}/data/clinicalAnnotation"
        params = {"location.genes.symbol": gene_symbol}
        retryable = {429, 500, 502, 503, 504}

        for attempt in range(2):
            try:
                timeout = aiohttp.ClientTimeout(total=60)
                async with aiohttp.ClientSession(timeout=timeout) as session:
                    async with session.get(search_url, params=params) as response:
                        if response.status == 200:
                            data = await response.json()
                            if not isinstance(data, dict) or not isinstance(data.get("data", []), list):
                                raise ValueError("ClinPGx returned an invalid response shape")
                            self._pharmgkb_last_status = "available"
                            for annotation in data.get("data", []):
                                loe = annotation.get("levelOfEvidence") or {}
                                level_term = loe.get("term")
                                variant_label = annotation.get("name")
                                acc = annotation.get("accessionId")
                                url = (f"https://www.pharmgkb.org/clinicalAnnotation/{acc}"
                                       if acc else None)
                                for chemical in annotation.get("relatedChemicals", []) or []:
                                    annotations.append({
                                        "gene": gene_symbol,
                                        "variant": variant_label,
                                        "drug": chemical.get("name"),
                                        "phenotype": annotation.get("phenotypeCategory"),
                                        "evidence_level": level_term,
                                        "score": annotation.get("score"),
                                        "url": url,
                                        "source": "pharmgkb",
                                    })
                            return annotations

                        if response.status not in retryable:
                            self._pharmgkb_last_status = "unavailable"
                            logger.warning("ClinPGx API error: %s", response.status)
                            return []
                        logger.warning("ClinPGx temporary API error: %s", response.status)
            except (aiohttp.ClientError, asyncio.TimeoutError, ValueError) as exc:
                logger.warning("ClinPGx request failed (attempt %s): %s", attempt + 1, exc)

            if attempt == 0:
                await asyncio.sleep(0.5)

        self._pharmgkb_last_status = "unavailable"
        return []
        
    def _combine_drug_results(self, cached_drugs: List[Dict[str, Any]], 
                             drugbank_drugs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Combine and deduplicate drug results from different sources"""
        
        # Use drug name as the key for deduplication
        combined_drugs = {}
        
        # Add cached drugs
        for drug in cached_drugs:
            key = drug.get('drug_name', '').lower()
            if key and key not in combined_drugs:
                combined_drugs[key] = drug
                
        # Add DrugBank drugs
        for drug in drugbank_drugs:
            key = drug.get('drug_name', '').lower()
            if key and key not in combined_drugs:
                combined_drugs[key] = drug
            elif key in combined_drugs:
                # Merge information from multiple sources
                existing = combined_drugs[key]
                existing['source'] = f"{existing.get('source', '')},drugbank_api".strip(',')
                if not existing.get('drug_id') and drug.get('drug_id'):
                    existing['drug_id'] = drug['drug_id']
                    
        return list(combined_drugs.values())
        
    def _rank_drug_matches(self, drugs: List[Dict[str, Any]], 
                          pharmgkb_annotations: List[Dict[str, Any]], 
                          hgvsp: Optional[str],
                          evidence_available: bool = True) -> List[Dict[str, Any]]:
        """Rank drug matches based on evidence and relevance"""
        
        # Create mapping of drug names to PharmGKB annotations
        pharmgkb_map = {}
        for annotation in pharmgkb_annotations:
            drug_name = annotation.get('drug', '').lower()
            if drug_name:
                if drug_name not in pharmgkb_map:
                    pharmgkb_map[drug_name] = []
                pharmgkb_map[drug_name].append(annotation)
                
        # Score and rank drugs
        scored_drugs = []

        # PharmGKB level -> (score, normalized label). Levels: 1A,1B,2A,2B,3,4
        level_scores = {
            "1a": (50, "1A"), "1b": (45, "1B"),
            "2a": (35, "2A"), "2b": (30, "2B"),
            "3": (20, "3"), "4": (10, "4"),
        }

        for drug in drugs:
            drug_name = drug.get('drug_name', '').lower()
            score = 10  # base score for having a drug-target relationship
            evidence_level = "no_pgx_evidence" if evidence_available else "evidence_unavailable"
            best_rank = -1

            if drug_name in pharmgkb_map:
                pgx_annotations = pharmgkb_map[drug_name]

                for annotation in pgx_annotations:
                    term = str(annotation.get('evidence_level') or '').lower().strip()
                    if term in level_scores:
                        add_score, label = level_scores[term]
                        if add_score > best_rank:
                            best_rank = add_score
                            score = 10 + add_score
                            evidence_level = f"PharmGKB Level {label}"
                    elif best_rank < 0:
                        # Annotation exists but level is unusual/unspecified
                        evidence_level = "PharmGKB (level n/a)"

                # Attach the PharmGKB annotations (with confirmation URLs) to the drug
                drug['pharmgkb_annotations'] = pgx_annotations
                # Surface the best confirmation URL directly on the drug
                for annotation in pgx_annotations:
                    if annotation.get('url'):
                        drug['pharmgkb_url'] = annotation['url']
                        break

            drug['evidence_score'] = score
            drug['evidence_level'] = evidence_level
            scored_drugs.append(drug)

        # Sort by score (descending)
        scored_drugs.sort(key=lambda x: x['evidence_score'], reverse=True)

        return scored_drugs
        
    async def cache_drug_data(self, drug_data: List[Dict[str, Any]]):
        """Cache drug-target relationships in local database"""
        
        try:
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()
            
            for drug in drug_data:
                cursor.execute('''
                    INSERT OR REPLACE INTO drug_targets 
                    (drug_id, drug_name, target_gene_symbol, target_uniprot_id, 
                     action, organism, pharmacology)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                ''', (
                    drug.get('drug_id'),
                    drug.get('drug_name'),
                    drug.get('target_gene_symbol'),
                    drug.get('target_uniprot_id'),
                    drug.get('action'),
                    drug.get('organism', 'Human'),
                    drug.get('pharmacology')
                ))
                
            conn.commit()
            conn.close()
            logger.info(f"Cached {len(drug_data)} drug records")
            
        except Exception as e:
            logger.error(f"Error caching drug data: {str(e)}")
            
    async def populate_known_pharmacogenes(self):
        """Populate database with known pharmacogenes and their drug interactions"""
        
        # Known pharmacogenes with key drug interactions
        known_pgx_genes = {
            'CYP2D6': [
                {'drug_name': 'Codeine', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Tramadol', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Metoprolol', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Risperidone', 'action': 'substrate', 'pharmacology': 'Metabolism'}
            ],
            'CYP2C19': [
                {'drug_name': 'Clopidogrel', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Omeprazole', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Escitalopram', 'action': 'substrate', 'pharmacology': 'Metabolism'}
            ],
            'CYP2C9': [
                {'drug_name': 'Warfarin', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Phenytoin', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Celecoxib', 'action': 'substrate', 'pharmacology': 'Metabolism'}
            ],
            'TPMT': [
                {'drug_name': 'Azathioprine', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Mercaptopurine', 'action': 'substrate', 'pharmacology': 'Metabolism'}
            ],
            'DPYD': [
                {'drug_name': 'Fluorouracil', 'action': 'substrate', 'pharmacology': 'Metabolism'},
                {'drug_name': 'Capecitabine', 'action': 'substrate', 'pharmacology': 'Metabolism'}
            ],
            'VKORC1': [
                {'drug_name': 'Warfarin', 'action': 'target', 'pharmacology': 'Vitamin K recycling'}
            ]
        }
        
        # Populate database
        drug_data = []
        for gene_symbol, drugs in known_pgx_genes.items():
            for drug in drugs:
                drug_data.append({
                    'drug_id': f"pgx_{drug['drug_name'].lower().replace(' ', '_')}",
                    'drug_name': drug['drug_name'],
                    'target_gene_symbol': gene_symbol,
                    'target_uniprot_id': None,
                    'action': drug['action'],
                    'organism': 'Human',
                    'pharmacology': drug['pharmacology']
                })
                
        await self.cache_drug_data(drug_data)
        logger.info("Populated known pharmacogenes database")