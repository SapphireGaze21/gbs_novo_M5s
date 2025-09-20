# market_analysis.py

# ==============================================================================
# Part 1: Data Acquisition
# Fetches scientific article data from the PubMed API.
# ==============================================================================
import pandas as pd
from Bio import Entrez

print("Starting Part 1: Data Acquisition from PubMed...")

# Always provide your email to the NCBI APIs
Entrez.email = "your.email@example.com" 

# Define the search query for PubMed
search_query = "obesity India prevalence"

# Use Entrez esearch to get a list of article IDs (PMIDs) that match the query
# retmax="100" limits the results to the first 100 articles
handle = Entrez.esearch(db="pubmed", term=search_query, retmax="100")
record = Entrez.read(handle)
handle.close()
id_list = record["IdList"]

# Use Entrez efetch to retrieve the full records for those IDs in XML format
handle = Entrez.efetch(db="pubmed", id=id_list, rettype="xml", retmode="text")
records = Entrez.read(handle)
handle.close()

# Parse the XML to extract the title and abstract for each article
titles = []
abstracts = []

if 'PubmedArticle' in records:
    for pubmed_article in records['PubmedArticle']:
        # Append title, with a check in case it's missing
        if 'ArticleTitle' in pubmed_article['MedlineCitation']['Article']:
            titles.append(pubmed_article['MedlineCitation']['Article']['ArticleTitle'])
        else:
            titles.append("No Title Found")
            
        # Append abstract, with a check in case it's missing
        if 'Abstract' in pubmed_article['MedlineCitation']['Article']:
            abstract_text = " ".join(pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText'])
            abstracts.append(abstract_text)
        else:
            abstracts.append("") # Use empty string if no abstract

# Store the raw data in a Pandas DataFrame
df_raw = pd.DataFrame({
    'title': titles,
    'abstract': abstracts
})

print(f"Successfully fetched {len(df_raw)} articles.")
print("Raw Data Head:")
print(df_raw.head())
print("-" * 50)


# ==============================================================================
# Part 2: Data Cleaning
# Uses the spaCy library to clean the fetched abstracts.
# ==============================================================================
import spacy
import re

print("\nStarting Part 2: Data Cleaning with spaCy...")

# Load the small English model for spaCy
# Make sure you have run: python -m spacy download en_core_web_sm
nlp = spacy.load("en_core_web_sm")

def clean_text(text):
    """
    Applies a basic NLP pipeline to clean text:
    - Converts to lowercase
    - Removes stop words and punctuation
    - Lemmatizes tokens (reduces words to their root form)
    """
    if not text:
        return ""
        
    doc = nlp(text.lower())
    
    # Create a list of lemmatized tokens that are not stop words or punctuation
    clean_tokens = []
    for token in doc:
        if not token.is_stop and not token.is_punct and token.is_alpha:
            clean_tokens.append(token.lemma_)
            
    return " ".join(clean_tokens)

# Apply the cleaning function to the 'abstract' column of our DataFrame
df_raw['cleaned_abstract'] = df_raw['abstract'].apply(clean_text)

print("Cleaning process complete.")
print("Final DataFrame Head with Cleaned Text:")
print(df_raw[['title', 'cleaned_abstract']].head())

# Optionally, save the final DataFrame to a CSV file
df_raw.to_csv('cleaned_pubmed_data.csv', index=False)
print("\nFinal data saved to 'cleaned_pubmed_data.csv'")