# AI-Powered Market Analysis for Wegovy in India 

This project, developed by **Team 'the m5s'**, presents a data-driven solution to inform the commercial strategy for Wegovy in India.We use an AI-powered platform to transform diverse, unstructured data from medical papers and national health surveys into an interactive analytics dashboard.

The platform provides crucial insights into obesity prevalence, patient profiles, and treatment patterns, enabling precise market segmentation, optimized pricing, and accurate demand forecasting for Novo Nordisk.


## Table of Contents

-   [Project Objective](#-project-objective)
-   [Key Features](#-key-features)
-   [Live Demo & Screenshots](#-live-demo--screenshots)
-   [Technical Architecture](#-technical-architecture)
-   [Key Insights & Commercial Strategy](#-key-insights--commercial-strategy)
-   [Getting Started](#-getting-started)
-   [Project Structure](#-project-structure)

## Project Objective

To conduct a comprehensive market analysis that quantifies obesity prevalence, patient profiles, and treatment patterns in India, providing data-driven insights to inform the commercial strategy for Wegovy.Our goal is to empower both Novo Nordisk with a clear go-to-market strategy and the community with data for targeted public health initiatives.

##  Key Features

Our Health Data Analytics Platform is built on three core pillars:

1. **ML Insights**: Delivers obesity prevalence forecasts with machine learning predictions and interactive visualizations across Indian states and districts.
2. **EDA Dashboard**: An exploratory data analysis tool for NFHS-5 health indicators, featuring interactive filtering and dynamic charts to explore demographic and health data.
3. **NLP Insights**: Leverages semantic search and AI-powered document summarization to analyze medical literature and extract insights that keyword searches would miss.


## Technical Architecture

We implemented a three-step pipeline to transform scattered information into a powerful strategic tool.

#### 1. Ingest & Clean
We gathered unstructured data from sources like **PubMed, NFHS, and WHO reports**. An NLP pipeline then cleaned and structured this information for analysis.
* **Tech Stack**: Python, Pandas, SpaCy (for NER and text cleaning).

#### 2. Analyze & Understand
We used advanced AI to cluster documents and perform semantic searches, uncovering hidden themes and relationships in the data.
* **Tech Stack**: Hugging Face (for embeddings), MongoDB (vector database), Scikit-learn (K-Means, DBSCAN).

#### 3. Predict & Visualize
We trained machine learning models to forecast future obesity trends and presented the insights in a dynamic, interactive dashboard.
* **Tech Stack**: Scikit-learn & TensorFlow (predictive modeling); React & D3.js (interactive dashboard).

## 💡 Key Insights & Commercial Strategy

Our analysis uncovered several key findings that directly inform a targeted commercial strategy for Wegovy.

### Key Insights

* **Geographic Hotspots**: Punjab, Delhi, and Chandigarh are predicted to have the highest obesity rates, exceeding 38%.
* **Urban-Rural Divide**: Predicted obesity among females is significantly higher in urban areas (33.8%) compared to rural areas (23.9%).
* **High Co-morbidity**: A strong correlation exists between high blood pressure and high blood sugar, identifying a clear patient profile for Wegovy.
* **Market Access Gaps**: A digital divide for women and low health insurance coverage are significant barriers to access.

### Recommended Commercial Strategy

1.  **Prioritize High-Potential Markets**: Focus sales and marketing efforts on top-tier urban hotspots like Delhi and Punjab.
2.  **Position for Co-morbid Patients**: Target healthcare professional education on Wegovy's benefits for patients with hypertension and other related risks.
3. **Optimize Pricing & Access**: Develop tiered pricing and Patient Assistance Programs (PAPs) to overcome affordability barriers.
4.  **Launch a Female-Focused Digital Program**: Create a mobile-first (WhatsApp/SMS) support program to bypass the digital access gap for women.
5.  **Execute a Phased Pilot Launch**: Engage Key Opinion Leaders (KOLs) in 3-5 high-potential states to build advocacy and guide a targeted rollout.

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
