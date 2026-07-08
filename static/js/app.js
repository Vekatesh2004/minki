// Pharmacogenomics Pipeline Web Interface - JavaScript

// Global variables
let systemStatus = 'unknown';
let analysisInProgress = false;

// Initialize when document is ready
$(document).ready(function() {
    initializeApp();
});

function initializeApp() {
    // Check system status
    checkSystemStatus();
    
    // Initialize tooltips
    if (typeof bootstrap !== 'undefined') {
        var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
        var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
            return new bootstrap.Tooltip(tooltipTriggerEl);
        });
    }
    
    // Initialize file upload handlers
    initializeFileUpload();
    
    // Set up periodic status checks
    setInterval(checkSystemStatus, 30000); // Check every 30 seconds
}

function checkSystemStatus() {
    $.ajax({
        url: '/api/status',
        method: 'GET',
        timeout: 5000,
        success: function(data) {
            systemStatus = data.status;
            updateStatusIndicator(data);
        },
        error: function() {
            systemStatus = 'error';
            updateStatusIndicator({status: 'error', pipeline_initialized: false});
        }
    });
}

function updateStatusIndicator(data) {
    const statusElement = $('#status-indicator');
    
    if (!statusElement.length) return;
    
    let statusHTML = '';
    let statusClass = '';
    
    if (data.status === 'active' && data.pipeline_initialized) {
        statusHTML = '<i class="fas fa-check-circle"></i> System Ready';
        statusClass = 'text-success';
    } else if (data.status === 'active') {
        statusHTML = '<i class="fas fa-clock"></i> Initializing...';
        statusClass = 'text-warning';
    } else {
        statusHTML = '<i class="fas fa-exclamation-triangle"></i> System Error';
        statusClass = 'text-danger';
    }
    
    statusElement.html(`<span class="${statusClass}">${statusHTML}</span>`);
}

function initializeFileUpload() {
    // File drag and drop support
    const fileInput = document.getElementById('vcf_file');
    
    if (fileInput) {
        // Prevent default drag behaviors
        ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
            fileInput.addEventListener(eventName, preventDefaults, false);
        });
        
        function preventDefaults(e) {
            e.preventDefault();
            e.stopPropagation();
        }
        
        // Highlight drop area when item is dragged over
        ['dragenter', 'dragover'].forEach(eventName => {
            fileInput.addEventListener(eventName, highlight, false);
        });
        
        ['dragleave', 'drop'].forEach(eventName => {
            fileInput.addEventListener(eventName, unhighlight, false);
        });
        
        function highlight(e) {
            fileInput.classList.add('border-primary');
        }
        
        function unhighlight(e) {
            fileInput.classList.remove('border-primary');
        }
        
        // Handle dropped files
        fileInput.addEventListener('drop', handleDrop, false);
        
        function handleDrop(e) {
            const dt = e.dataTransfer;
            const files = dt.files;
            
            if (files.length > 0) {
                fileInput.files = files;
                validateFile(files[0]);
            }
        }
    }
}

function validateFile(file) {
    const maxSize = 50 * 1024 * 1024; // 50MB
    const allowedExtensions = ['vcf', 'gz'];
    
    // Check file size
    if (file.size > maxSize) {
        showAlert('File size exceeds 50MB limit', 'danger');
        return false;
    }
    
    // Check file extension
    const fileName = file.name.toLowerCase();
    const hasValidExtension = allowedExtensions.some(ext => 
        fileName.endsWith('.' + ext) || fileName.endsWith('.vcf.gz')
    );
    
    if (!hasValidExtension) {
        showAlert('Please select a VCF file (.vcf or .vcf.gz)', 'danger');
        return false;
    }
    
    return true;
}

function showAlert(message, type = 'info', timeout = 5000) {
    const alertHTML = `
        <div class="alert alert-${type} alert-dismissible fade show" role="alert">
            <i class="fas fa-${getAlertIcon(type)}"></i>
            ${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        </div>
    `;
    
    // Find or create alert container
    let alertContainer = $('#alert-container');
    if (!alertContainer.length) {
        $('main.container').prepend('<div id="alert-container"></div>');
        alertContainer = $('#alert-container');
    }
    
    alertContainer.prepend(alertHTML);
    
    // Auto-remove after timeout
    if (timeout > 0) {
        setTimeout(() => {
            alertContainer.find('.alert').first().alert('close');
        }, timeout);
    }
}

function getAlertIcon(type) {
    const icons = {
        'success': 'check-circle',
        'danger': 'exclamation-triangle',
        'warning': 'exclamation-triangle',
        'info': 'info-circle'
    };
    return icons[type] || 'info-circle';
}

// Analysis utilities
function startAnalysis(filename, sampleId, analysisType = 'basic') {
    if (analysisInProgress) {
        showAlert('Analysis already in progress', 'warning');
        return;
    }
    
    analysisInProgress = true;
    
    return $.ajax({
        url: '/api/analyze',
        method: 'POST',
        contentType: 'application/json',
        data: JSON.stringify({
            filename: filename,
            sample_id: sampleId,
            analysis_type: analysisType
        }),
        timeout: 300000, // 5 minute timeout
        success: function(response) {
            analysisInProgress = false;
            return response;
        },
        error: function(xhr, status, error) {
            analysisInProgress = false;
            throw error;
        }
    });
}

// Utility functions
function formatFileSize(bytes) {
    if (bytes === 0) return '0 Bytes';
    
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}

function formatPercentage(value) {
    return (value * 100).toFixed(1) + '%';
}

function formatNumber(num) {
    return num.toLocaleString();
}

function copyToClipboard(text) {
    if (navigator.clipboard) {
        navigator.clipboard.writeText(text).then(() => {
            showAlert('Copied to clipboard', 'success', 2000);
        });
    } else {
        // Fallback for older browsers
        const textArea = document.createElement('textarea');
        textArea.value = text;
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand('copy');
        document.body.removeChild(textArea);
        showAlert('Copied to clipboard', 'success', 2000);
    }
}

// Download utilities
function downloadJSON(data, filename) {
    const dataStr = JSON.stringify(data, null, 2);
    const dataUri = 'data:application/json;charset=utf-8,'+ encodeURIComponent(dataStr);
    
    const exportFileDefaultName = filename || 'data.json';
    
    const linkElement = document.createElement('a');
    linkElement.setAttribute('href', dataUri);
    linkElement.setAttribute('download', exportFileDefaultName);
    linkElement.click();
}

function downloadCSV(data, filename) {
    let csvContent = '';
    
    // Assume data is array of objects
    if (data.length > 0) {
        // Header
        const headers = Object.keys(data[0]);
        csvContent += headers.join(',') + '\n';
        
        // Rows
        data.forEach(row => {
            const values = headers.map(header => {
                const value = row[header];
                return typeof value === 'string' ? `"${value.replace(/"/g, '""')}"` : value;
            });
            csvContent += values.join(',') + '\n';
        });
    }
    
    const dataUri = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csvContent);
    const exportFileDefaultName = filename || 'data.csv';
    
    const linkElement = document.createElement('a');
    linkElement.setAttribute('href', dataUri);
    linkElement.setAttribute('download', exportFileDefaultName);
    linkElement.click();
}

// Progress tracking
function createProgressTracker(steps) {
    return {
        currentStep: 0,
        steps: steps,
        
        nextStep: function() {
            if (this.currentStep < this.steps.length - 1) {
                this.currentStep++;
            }
            return this.getCurrentStep();
        },
        
        getCurrentStep: function() {
            return this.steps[this.currentStep];
        },
        
        getProgress: function() {
            return ((this.currentStep + 1) / this.steps.length) * 100;
        },
        
        isComplete: function() {
            return this.currentStep === this.steps.length - 1;
        }
    };
}

// Table utilities
function createDataTable(data, containerId, options = {}) {
    const container = $(`#${containerId}`);
    
    if (!data || data.length === 0) {
        container.html('<p class="text-muted">No data available</p>');
        return;
    }
    
    const headers = Object.keys(data[0]);
    
    let tableHTML = `
        <div class="table-responsive">
            <table class="table table-striped table-hover">
                <thead class="table-dark">
                    <tr>
                        ${headers.map(header => `<th>${header}</th>`).join('')}
                    </tr>
                </thead>
                <tbody>
    `;
    
    const maxRows = options.maxRows || data.length;
    data.slice(0, maxRows).forEach(row => {
        tableHTML += '<tr>';
        headers.forEach(header => {
            const value = row[header];
            const displayValue = options.formatters && options.formatters[header] 
                ? options.formatters[header](value) 
                : (value || 'N/A');
            tableHTML += `<td>${displayValue}</td>`;
        });
        tableHTML += '</tr>';
    });
    
    tableHTML += `
                </tbody>
            </table>
        </div>
    `;
    
    if (data.length > maxRows) {
        tableHTML += `<p class="text-muted">Showing ${maxRows} of ${data.length} rows</p>`;
    }
    
    container.html(tableHTML);
}

// Chart utilities (if Chart.js is available)
function createChart(canvasId, type, data, options = {}) {
    if (typeof Chart === 'undefined') {
        console.warn('Chart.js not loaded');
        return;
    }
    
    const ctx = document.getElementById(canvasId);
    if (!ctx) {
        console.error(`Canvas element ${canvasId} not found`);
        return;
    }
    
    return new Chart(ctx, {
        type: type,
        data: data,
        options: {
            responsive: true,
            maintainAspectRatio: false,
            ...options
        }
    });
}

// Error handling
window.addEventListener('error', function(e) {
    console.error('JavaScript error:', e.error);
    // Don't show alert for every JS error to avoid spam
});

// Handle AJAX errors globally
$(document).ajaxError(function(event, xhr, settings, error) {
    if (xhr.status === 0) {
        showAlert('Network error - please check your connection', 'danger');
    } else if (xhr.status === 500) {
        showAlert('Server error - please try again later', 'danger');
    } else if (xhr.status === 404) {
        showAlert('Resource not found', 'danger');
    }
});

// Export functions for global use
window.PharmacogenomicsApp = {
    showAlert,
    startAnalysis,
    downloadJSON,
    downloadCSV,
    formatFileSize,
    formatPercentage,
    formatNumber,
    copyToClipboard,
    createDataTable,
    createProgressTracker
};