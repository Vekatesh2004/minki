# Deploying the Pharmacogenomics Pipeline to AWS EC2

This guide deploys `simple_backend.py` (the working FastAPI app) on a single
Ubuntu EC2 instance, kept alive by systemd and served through nginx on port 80.

Architecture:
```
Browser  ->  EC2 :80 (nginx)  ->  127.0.0.1:8000 (uvicorn / FastAPI)  ->  Ensembl VEP / UniProt / AlphaFold / PharmGKB
```

---

## 1. Launch an EC2 instance

1. AWS Console -> EC2 -> **Launch instance**.
2. **Name**: pharmacogenomics
3. **AMI**: Ubuntu Server 22.04 LTS (64-bit x86).
4. **Instance type**: `t3.small` (2 GB RAM) minimum; `t3.medium` (4 GB) recommended.
   The app calls external APIs and parses large VCFs, so avoid `t2.micro`.
5. **Key pair**: create or select one (e.g. `minki-key.pem`). Download and keep it safe.
6. **Network / Security group** — add inbound rules:
   - SSH (TCP 22) — Source: **My IP** (not 0.0.0.0/0)
   - HTTP (TCP 80) — Source: Anywhere (0.0.0.0/0)
   - (Add HTTPS TCP 443 later if you set up a domain + TLS)
7. **Storage**: 20 GB gp3 is plenty.
8. Launch.

Note the instance's **Public IPv4 address** (e.g. `13.234.56.78`).

---

## 2. Connect via SSH

From your local machine, in the folder with your key:

```bash
chmod 400 minki-key.pem
ssh -i minki-key.pem ubuntu@YOUR_PUBLIC_IP
```

---

## 3. Install system packages

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y python3-venv python3-pip nginx git
```

---

## 4. Get the code onto the server

**Option A — git (if your project is in a repo):**
```bash
cd ~
git clone YOUR_REPO_URL minki
```

**Option B — copy from your laptop with scp** (run this on your LAPTOP, not the server):
```bash
# from /home/venkatesh-g/Documents
scp -i minki-key.pem -r minki ubuntu@YOUR_PUBLIC_IP:/home/ubuntu/minki
```

---

## 5. Set up the Python environment

On the server:

```bash
cd ~/minki
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements_deploy.txt
```

Quick smoke test (Ctrl+C after you see "Uvicorn running"):
```bash
HOST=127.0.0.1 PORT=8000 RELOAD=0 python simple_backend.py
```

---

## 6. Run it as a systemd service (auto-start + restart)

```bash
sudo cp ~/minki/deploy/pharmacogenomics.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable pharmacogenomics
sudo systemctl start pharmacogenomics

# check it's running
sudo systemctl status pharmacogenomics
curl -s http://127.0.0.1:8000/health
```

If `status` shows errors, view logs with:
```bash
sudo journalctl -u pharmacogenomics -f
```

---

## 7. Put nginx in front (port 80)

```bash
sudo cp ~/minki/deploy/nginx-pharmacogenomics.conf /etc/nginx/sites-available/pharmacogenomics
sudo ln -s /etc/nginx/sites-available/pharmacogenomics /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default   # remove the default welcome page
sudo nginx -t                                 # test config
sudo systemctl restart nginx
```

Now open in your browser:
```
http://YOUR_PUBLIC_IP
```
You should see the Pharmacogenomics Pipeline UI. Upload
`examples/sample_pharmacogenomics.vcf` to test.

---

## 8. (Optional) Domain + HTTPS

If you have a domain, point an A record at the EC2 public IP, then:
```bash
sudo apt install -y certbot python3-certbot-nginx
sudo certbot --nginx -d yourdomain.com
```
Certbot edits the nginx config and auto-renews the certificate.

---

## Updating the app later

```bash
ssh -i minki-key.pem ubuntu@YOUR_PUBLIC_IP
cd ~/minki
git pull            # or scp the changed files again
source venv/bin/activate
pip install -r requirements_deploy.txt   # only if deps changed
sudo systemctl restart pharmacogenomics
```

---

## Tuning notes

- **Annotate more variants**: edit `MAX_VARIANTS_TO_VEP` in
  `/etc/systemd/system/pharmacogenomics.service`, then
  `sudo systemctl daemon-reload && sudo systemctl restart pharmacogenomics`.
  Higher values = more Ensembl VEP API calls = slower.
- **Costs**: a `t3.small` runs ~US$15/month if left on 24/7. Stop the instance
  when not in use to save money (the public IP changes on stop/start unless you
  attach an Elastic IP).

---

## Security reminders

- Keep SSH (port 22) restricted to your IP.
- This app has **no authentication** — anyone who knows the IP can use it. For a
  private tool, restrict port 80 in the security group to your IP, or add HTTP
  basic auth in nginx.
- Don't commit real API keys into `config.json` in a public repo.
