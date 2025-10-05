# AI-Powered Market Analysis for Wegovy in India 

[cite_start]This project, developed by **Team 'the m5s'**, presents a data-driven solution to inform the commercial strategy for Wegovy in India[cite: 1, 2]. [cite_start]We use an AI-powered platform to transform diverse, unstructured data from medical papers and national health surveys into an interactive analytics dashboard[cite: 11, 12].

[cite_start]The platform provides crucial insights into obesity prevalence, patient profiles, and treatment patterns, enabling precise market segmentation, optimized pricing, and accurate demand forecasting for Novo Nordisk[cite: 12, 17].


## Table of Contents

-   [Project Objective](#-project-objective)
-   [Key Features](#-key-features)
-   [Live Demo & Screenshots](#-live-demo--screenshots)
-   [Technical Architecture](#-technical-architecture)
-   [Key Insights & Commercial Strategy](#-key-insights--commercial-strategy)
-   [Getting Started](#-getting-started)
-   [Project Structure](#-project-structure)

## Project Objective

[cite_start]To conduct a comprehensive market analysis that quantifies obesity prevalence, patient profiles, and treatment patterns in India, providing data-driven insights to inform the commercial strategy for Wegovy[cite: 1]. [cite_start]Our goal is to empower both Novo Nordisk with a clear go-to-market strategy and the community with data for targeted public health initiatives[cite: 17, 18].

##  Key Features

Our Health Data Analytics Platform is built on three core pillars:

1.  [cite_start]**🧠 ML Insights**: Delivers obesity prevalence forecasts with machine learning predictions and interactive visualizations across Indian states and districts[cite: 29, 56].
2.  [cite_start]**📊 EDA Dashboard**: An exploratory data analysis tool for NFHS-5 health indicators, featuring interactive filtering and dynamic charts to explore demographic and health data[cite: 56].
3.  [cite_start]**🔎 NLP Insights**: Leverages semantic search and AI-powered document summarization to analyze medical literature and extract insights that keyword searches would miss[cite: 26, 56].

## Live Demo & Screenshots

Here is a preview of the interactive dashboard, which serves as the primary interface for exploring the project's findings.

[cite_start]*The main landing page of the Health Data Analytics Platform.* [cite: 56]

---

[cite_start]*ML-powered forecasts showing predicted obesity prevalence across 107 districts.* [cite: 56]

---

[cite_start]*Interactive charts in the EDA dashboard exploring the prevalence of obesity by state.* [cite: 56]

## 🛠️ Technical Architecture

[cite_start]We implemented a three-step pipeline to transform scattered information into a powerful strategic tool[cite: 20].

#### 1. Ingest & Clean
[cite_start]We gathered unstructured data from sources like **PubMed, NFHS, and WHO reports**[cite: 22]. [cite_start]An NLP pipeline then cleaned and structured this information for analysis[cite: 23].
* [cite_start]**Tech Stack**: Python, Pandas, SpaCy (for NER and text cleaning)[cite: 24].

#### 2. Analyze & Understand
[cite_start]We used advanced AI to cluster documents and perform semantic searches, uncovering hidden themes and relationships in the data[cite: 26].
* [cite_start]**Tech Stack**: Hugging Face (for embeddings), MongoDB (vector database), Scikit-learn (K-Means, DBSCAN)[cite: 27].

#### 3. Predict & Visualize
[cite_start]We trained machine learning models to forecast future obesity trends and presented the insights in a dynamic, interactive dashboard[cite: 29].
* **Tech Stack**: Scikit-learn & TensorFlow (predictive modeling); [cite_start]React & D3.js (interactive dashboard)[cite: 30].

## 💡 Key Insights & Commercial Strategy

Our analysis uncovered several key findings that directly inform a targeted commercial strategy for Wegovy.

### Key Insights

* [cite_start]**Geographic Hotspots**: Punjab, Delhi, and Chandigarh are predicted to have the highest obesity rates, exceeding 38%[cite: 35].
* [cite_start]**Urban-Rural Divide**: Predicted obesity among females is significantly higher in urban areas (33.8%) compared to rural areas (23.9%)[cite: 36].
* [cite_start]**High Co-morbidity**: A strong correlation exists between high blood pressure and high blood sugar, identifying a clear patient profile for Wegovy[cite: 38, 39].
* [cite_start]**Market Access Gaps**: A digital divide for women and low health insurance coverage are significant barriers to access[cite: 41, 42].

### Recommended Commercial Strategy

1.  [cite_start]**Prioritize High-Potential Markets**: Focus sales and marketing efforts on top-tier urban hotspots like Delhi and Punjab[cite: 45].
2.  [cite_start]**Position for Co-morbid Patients**: Target healthcare professional education on Wegovy's benefits for patients with hypertension and other related risks[cite: 47].
3.  [cite_start]**Optimize Pricing & Access**: Develop tiered pricing and Patient Assistance Programs (PAPs) to overcome affordability barriers[cite: 49].
4.  [cite_start]**Launch a Female-Focused Digital Program**: Create a mobile-first (WhatsApp/SMS) support program to bypass the digital access gap for women[cite: 51].
5.  [cite_start]**Execute a Phased Pilot Launch**: Engage Key Opinion Leaders (KOLs) in 3-5 high-potential states to build advocacy and guide a targeted rollout[cite: 53].

## Getting Started

To get a local copy up and running, follow these simple steps.

### Prerequisites

* Node.js & npm (for frontend)
* Python 3.8+ & pip (for backend)
* MongoDB instance

### Installation

1.  **Clone the repo**
    ```sh
    git clone [https://github.com/SapphireGaze21/gbs_novo_M5s.git](https://github.com/SapphireGaze21/gbs_novo_M5s.git)
    cd gbs_novo_M5s/novo
    ```
2.  **Install Frontend Dependencies**
    ```sh
    cd frontend
    npm install
    ```
3.  **Install Backend Dependencies**
    ```sh
    cd ../backend
    pip install -r requirements.txt
    ```
4.  **Set up Environment Variables**
    Create a `.env` file in the `backend` directory and add your MongoDB connection string and any other necessary API keys.
    ```
    MONGO_URI="your_connection_string"
    ```
5.  **Run the application**
    * Start the backend server: `python app.py`
    * Start the frontend development server: `npm start`

## Project Structure

```
novo/
├── app/              # Main application routes, pages, and layouts
├── components/       # Reusable UI components (e.g., charts, cards)
├── lib/              # Utility functions and helper scripts
├── .gitignore        # Files and folders to be ignored by Git
├── next.config.mjs   # Next.js configuration file
├── package.json      # Project dependencies and scripts
├── tailwind.config.ts# Configuration for Tailwind CSS
└── README.md         # This file
```
