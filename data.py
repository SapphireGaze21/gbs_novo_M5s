# ==============================================================================
# market_analysis.py
# Robust PubMed Data Acquisition & Cleaning Script with Retry & Batching
# ==============================================================================

import pandas as pd
from Bio import Entrez
import spacy
import time
from http.client import IncompleteRead

print("Step 1: Searching PubMed and Fetching Articles")

# Always provide your email to NCBI
Entrez.email = "your.email@example.com"

# Set your search keywords here
search_query = "obesity India prevalence"

# Search PubMed articles - get up to 100 PMIDs
search_handle = Entrez.esearch(db="pubmed", term=search_query, retmax=100)
search_record = Entrez.read(search_handle)
search_handle.close()
id_list = search_record["IdList"]

print(f"Found {len(id_list)} articles matching query.")

# Function to fetch articles in batches with retry on failures
def fetch_pubmed_records(ids, batch_size=20, max_retries=3):
    all_records = []
    for start in range(0, len(ids), batch_size):
        batch_ids = ids[start:start+batch_size]
        retry = 0
        while retry < max_retries:
            try:
                handle = Entrez.efetch(db="pubmed", id=batch_ids, rettype="xml", retmode="text")
                records = Entrez.read(handle)
                handle.close()
                all_records.extend(records.get('PubmedArticle', []))
                # Success, break out of retry loop
                break
            except IncompleteRead as e:
                retry += 1
                print(f"IncompleteRead error, retry {retry}/{max_retries}...")
                time.sleep(2)  # short delay before retry
            except Exception as e:
                print(f"Unexpected error: {e}")
                break
        else:
            print(f"Failed to fetch batch starting at index {start} after {max_retries} retries.")
    return all_records

# Fetch the article details robustly
pubmed_articles = fetch_pubmed_records(id_list, batch_size=20)

print(f"Total articles fetched: {len(pubmed_articles)}")

# Extract titles and abstracts
titles = []
abstracts = []
for article in pubmed_articles:
    art = article['MedlineCitation']['Article']
    titles.append(art.get('ArticleTitle', 'No Title Found'))
    if 'Abstract' in art:
        abstracts.append(" ".join(art['Abstract']['AbstractText']))
    else:
        abstracts.append("")

# Create DataFrame from extracted data
df_raw = pd.DataFrame({'title': titles, 'abstract': abstracts})
print("Sample raw data:")
print(df_raw.head())
print("-" * 50)

print("\nStep 2: Cleaning abstracts with spaCy")

# Load spaCy English model (install with: python -m spacy download en_core_web_sm)
nlp = spacy.load("en_core_web_sm")

# Function to clean text using spaCy
def clean_text(text):
    if not text:
        return ""
    doc = nlp(text.lower())
    clean_tokens = [token.lemma_ for token in doc if token.is_alpha and not token.is_stop and not token.is_punct]
    return " ".join(clean_tokens)

# Clean abstracts
df_raw['cleaned_abstract'] = df_raw['abstract'].apply(clean_text)

print("Sample cleaned data:")
print(df_raw[['title', 'cleaned_abstract']].head())

# Save cleaned data to CSV
df_raw.to_csv('cleaned_pubmed_data.csv', index=False)
print("Data saved to 'cleaned_pubmed_data.csv'.")
