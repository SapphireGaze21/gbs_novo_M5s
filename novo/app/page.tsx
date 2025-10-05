"use client"

import Link from "next/link"
import { Card } from "@/components/ui/card"
import { BarChart3, Search, Brain } from "lucide-react"

export default function HomePage() {
  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="container mx-auto px-6 py-12">
        <div className="text-center mb-12">
          <h1 className="text-4xl font-bold mb-4 bg-gradient-to-r from-cyan-400 to-teal-400 bg-clip-text text-transparent">
            Health Data Analytics Platform
          </h1>
          <p className="text-slate-300 text-lg max-w-2xl mx-auto">
            Explore comprehensive health insights through machine learning predictions, exploratory data analysis, and
            natural language processing.
          </p>
        </div>

        <div className="grid md:grid-cols-3 gap-8 max-w-6xl mx-auto">
          <Link href="/ml-insights">
            <Card className="p-8 bg-slate-900/50 backdrop-blur-sm border-slate-700/50 hover:border-cyan-400/50 transition-all duration-300 cursor-pointer group">
              <div className="text-center">
                <div className="w-16 h-16 mx-auto mb-4 bg-gradient-to-br from-cyan-400 to-teal-500 rounded-full flex items-center justify-center group-hover:scale-110 transition-transform">
                  <BarChart3 className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-semibold mb-2 text-white">ML Insights</h3>
                <p className="text-slate-400">
                  Obesity prevalence forecasting with machine learning predictions and interactive visualizations.
                </p>
              </div>
            </Card>
          </Link>

          <Link href="/eda">
            <Card className="p-8 bg-slate-900/50 backdrop-blur-sm border-slate-700/50 hover:border-cyan-400/50 transition-all duration-300 cursor-pointer group">
              <div className="text-center">
                <div className="w-16 h-16 mx-auto mb-4 bg-gradient-to-br from-cyan-400 to-teal-500 rounded-full flex items-center justify-center group-hover:scale-110 transition-transform">
                  <Search className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-semibold mb-2 text-white">EDA Dashboard</h3>
                <p className="text-slate-400">
                  Exploratory data analysis of NFHS-5 health indicators with interactive filtering and charts.
                </p>
              </div>
            </Card>
          </Link>

          <Link href="/nlp-insights">
            <Card className="p-8 bg-slate-900/50 backdrop-blur-sm border-slate-700/50 hover:border-cyan-400/50 transition-all duration-300 cursor-pointer group">
              <div className="text-center">
                <div className="w-16 h-16 mx-auto mb-4 bg-gradient-to-br from-cyan-400 to-teal-500 rounded-full flex items-center justify-center group-hover:scale-110 transition-transform">
                  <Brain className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-semibold mb-2 text-white">NLP Insights</h3>
                <p className="text-slate-400">
                  Medical literature analysis with semantic search and AI-powered document summarization.
                </p>
              </div>
            </Card>
          </Link>
        </div>
      </div>
    </div>
  )
}
