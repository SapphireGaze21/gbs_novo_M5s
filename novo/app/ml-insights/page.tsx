"use client"

import { useState, useEffect } from "react"
import { DashboardHeader } from "@/components/dashboard-header"
import { KpiCard } from "@/components/kpi-card"
import { GlassTable } from "@/components/glass-table"
import { GlassChart } from "@/components/glass-chart"
import { GlassCard } from "@/components/glass-card"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { Label } from "@/components/ui/label"
import { Slider } from "@/components/ui/slider"
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from "recharts"

interface PredictionData {
  state: string
  area: string
  pred_women_overweight_obese_pct: number
  pred_men_overweight_obese_pct: number
  [key: string]: any
}

export default function MLInsightsPage() {
  const [data, setData] = useState<PredictionData[]>([])
  const [filteredData, setFilteredData] = useState<PredictionData[]>([])
  const [selectedState, setSelectedState] = useState<string>("all")
  const [selectedArea, setSelectedArea] = useState<string>("all")
  const [prevalenceRange, setPrevalenceRange] = useState([0, 100])
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await fetch(
          "https://hebbkx1anhila5yf.public.blob.vercel-storage.com/all_predictions_multi-dDbWyfQC1C9nMCnrRvpMH2wePsitAl.csv",
        )
        const csvText = await response.text()

        const lines = csvText.split("\n")
        const headers = lines[0].split(",")

        const parsedData = lines
          .slice(1)
          .filter((line) => line.trim())
          .map((line) => {
            const values = line.split(",")
            const row: any = {}
            headers.forEach((header, index) => {
              const cleanHeader = header.trim()
              const value = values[index]?.trim()

              if (cleanHeader.includes("pred_") && !isNaN(Number(value))) {
                row[cleanHeader] = Number(value)
              } else {
                row[cleanHeader] = value
              }
            })
            return row
          })

        setData(parsedData)
        setFilteredData(parsedData)
        setLoading(false)
      } catch (error) {
        console.error("Error fetching data:", error)
        setLoading(false)
      }
    }

    fetchData()
  }, [])

  useEffect(() => {
    let filtered = data

    if (selectedState !== "all") {
      filtered = filtered.filter((row) => row.state === selectedState)
    }

    if (selectedArea !== "all") {
      filtered = filtered.filter((row) => row.area.toLowerCase() === selectedArea.toLowerCase())
    }

    filtered = filtered.filter((row) => {
      const womenPrev = row.pred_women_overweight_obese_pct || 0
      const menPrev = row.pred_men_overweight_obese_pct || 0
      const avgPrev = (womenPrev + menPrev) / 2
      return avgPrev >= prevalenceRange[0] && avgPrev <= prevalenceRange[1]
    })

    setFilteredData(filtered)
  }, [data, selectedState, selectedArea, prevalenceRange])

  const uniqueStates = [...new Set(data.map((row) => row.state))].sort()

  const kpiData = {
    totalDistricts: filteredData.length,
    avgPrevalence:
      filteredData.length > 0
        ? (
            filteredData.reduce(
              (sum, row) => sum + (row.pred_women_overweight_obese_pct + row.pred_men_overweight_obese_pct) / 2,
              0,
            ) / filteredData.length
          ).toFixed(1)
        : "0",
    highestRisk:
      filteredData.length > 0
        ? filteredData.reduce((max, row) => {
            const avgPrev = (row.pred_women_overweight_obese_pct + row.pred_men_overweight_obese_pct) / 2
            const maxAvg = (max.pred_women_overweight_obese_pct + max.pred_men_overweight_obese_pct) / 2
            return avgPrev > maxAvg ? row : max
          }).state
        : "N/A",
  }

  const chartData = filteredData
    .sort((a, b) => {
      const aAvg = (a.pred_women_overweight_obese_pct + a.pred_men_overweight_obese_pct) / 2
      const bAvg = (b.pred_women_overweight_obese_pct + b.pred_men_overweight_obese_pct) / 2
      return bAvg - aAvg
    })
    .slice(0, 15)
    .map((row) => ({
      name: `${row.state} (${row.area})`,
      prevalence: ((row.pred_women_overweight_obese_pct + row.pred_men_overweight_obese_pct) / 2).toFixed(1),
    }))

  const tableColumns = [
    { key: "state", label: "State/UT", sortable: true },
    { key: "area", label: "Area", sortable: true },
    {
      key: "pred_women_overweight_obese_pct",
      label: "Women Overweight/Obese (%)",
      sortable: true,
    },
    {
      key: "pred_men_overweight_obese_pct",
      label: "Men Overweight/Obese (%)",
      sortable: true,
    },
  ]

  const tableData = filteredData.slice(0, 50).map((row) => ({
    state: row.state,
    area: row.area,
    pred_women_overweight_obese_pct: row.pred_women_overweight_obese_pct?.toFixed(1) || "N/A",
    pred_men_overweight_obese_pct: row.pred_men_overweight_obese_pct?.toFixed(1) || "N/A",
  }))

  if (loading) {
    return (
      <div className="min-h-screen bg-slate-950 text-white p-6">
        <div className="container mx-auto">
          <div className="flex items-center justify-center h-64">
            <div className="text-cyan-400">Loading ML insights...</div>
          </div>
        </div>
      </div>
    )
  }

  return (
    <div className="min-h-screen bg-slate-950 text-white p-6">
      <div className="container mx-auto">
        <DashboardHeader
          title="ML Model Insights: Obesity Prevalence Forecast"
          subtitle="Machine learning predictions for overweight and obesity prevalence across Indian states and districts"
        />

        {/* KPI Cards */}
        <div className="grid md:grid-cols-3 gap-6 mb-8">
          <KpiCard title="Total Districts Analyzed" value={kpiData.totalDistricts} />
          <KpiCard title="Average Predicted Prevalence" value={`${kpiData.avgPrevalence}%`} />
          <KpiCard title="Highest Risk State" value={kpiData.highestRisk} />
        </div>

        {/* Filter Controls */}
        <GlassCard className="p-6 mb-8">
          <h3 className="text-lg font-semibold text-white mb-4">Filter Controls</h3>
          <div className="grid md:grid-cols-3 gap-6">
            <div>
              <Label className="text-slate-300 mb-2 block">State/UT</Label>
              <Select value={selectedState} onValueChange={setSelectedState}>
                <SelectTrigger className="bg-slate-800/50 border-slate-600 text-white">
                  <SelectValue />
                </SelectTrigger>
                <SelectContent className="bg-slate-800 border-slate-600">
                  <SelectItem value="all" className="text-white hover:bg-slate-700">
                    All States
                  </SelectItem>
                  {uniqueStates.map((state) => (
                    <SelectItem key={state} value={state} className="text-white hover:bg-slate-700">
                      {state}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>

            <div>
              <Label className="text-slate-300 mb-2 block">Area</Label>
              <Select value={selectedArea} onValueChange={setSelectedArea}>
                <SelectTrigger className="bg-slate-800/50 border-slate-600 text-white">
                  <SelectValue />
                </SelectTrigger>
                <SelectContent className="bg-slate-800 border-slate-600">
                  <SelectItem value="all" className="text-white hover:bg-slate-700">
                    All Areas
                  </SelectItem>
                  <SelectItem value="urban" className="text-white hover:bg-slate-700">
                    Urban
                  </SelectItem>
                  <SelectItem value="rural" className="text-white hover:bg-slate-700">
                    Rural
                  </SelectItem>
                </SelectContent>
              </Select>
            </div>

            <div>
              <Label className="text-slate-300 mb-2 block">
                Predicted Prevalence Range: {prevalenceRange[0]}% - {prevalenceRange[1]}%
              </Label>
              <Slider
                value={prevalenceRange}
                onValueChange={setPrevalenceRange}
                max={100}
                min={0}
                step={1}
                className="w-full"
              />
            </div>
          </div>
        </GlassCard>

        {/* Main Content Grid */}
        <div className="grid lg:grid-cols-2 gap-8">
          {/* Data Table */}
          <div>
            <GlassTable columns={tableColumns} data={tableData} className="h-fit" />
          </div>

          {/* Visualization Chart */}
          <div>
            <GlassChart title="Top 15 Districts by Predicted Prevalence">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={chartData} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                  <XAxis dataKey="name" stroke="#94a3b8" angle={-45} textAnchor="end" height={80} fontSize={10} />
                  <YAxis stroke="#94a3b8" />
                  <Tooltip
                    contentStyle={{
                      backgroundColor: "#1e293b",
                      border: "1px solid #06b6d4",
                      borderRadius: "8px",
                      color: "#fff",
                    }}
                  />
                  <Bar dataKey="prevalence" fill="#06b6d4" />
                </BarChart>
              </ResponsiveContainer>
            </GlassChart>
          </div>
        </div>
      </div>
    </div>
  )
}
