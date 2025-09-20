import pdfplumber
import os

pdf_files = [
    'stifeladaobesitydrugreview_07.01.2023.pdf',
    'epocrates-the-growing-glp-1-market-landscape-key-insights-for-pharma-marketers.pdf',
    'future-outlook-AOM-v3.pdf',
    'stifelobesityupdate_july2024.pdf',
    'Affordable-Access-to-GLP-1-Obesity-Medications-_-ICER-White-Paper-_-04.09.2025.pdf'
]

output_folder = 'extracted_texts'
os.makedirs(output_folder, exist_ok=True)

for pdf_file in pdf_files:
    full_text = ""
    with pdfplumber.open(pdf_file) as pdf:
        for page in pdf.pages:
            text = page.extract_text()
            if text:
                full_text += text + "\n\n"
    output_filepath = os.path.join(output_folder, pdf_file.replace('.pdf', '.txt'))
    with open(output_filepath, 'w', encoding='utf-8') as f:
        f.write(full_text)
