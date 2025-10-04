"use client"

import { useState, useEffect } from "react"
import { DashboardHeader } from "@/components/dashboard-header"
import { FilterSidebar } from "@/components/filter-sidebar"
import { GlassChart } from "@/components/glass-chart"
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, PieChart, Pie, Cell } from "recharts"

interface NFHSData {
  state: string
  area: string
  women_overweight_obese_pct: number
  men_overweight_obese_pct: number
  women_high_blood_pressure_pct: number
  men_high_blood_pressure_pct: number
  [key: string]: any
}

export default function EDAPage() {
  const [data, setData] = useState<NFHSData[]>([])
  const [filteredData, setFilteredData] = useState<NFHSData[]>([])
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    const fetchData = async () => {
      try {
        const response = await fetch(
          "https://hebbkx1anhila5yf.public.blob.vercel-storage.com/NFHS_5_cleaned_data-LmgKCIay6QqJ8YMdrsT5P2maEonf0b.csv",
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

              if (cleanHeader.includes("_pct") && !isNaN(Number(value))) {
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

  const handleFiltersChange = (filters: any) => {
    let filtered = data

    if (filters.states.length > 0) {
      filtered = filtered.filter((row) => filters.states.includes(row.state))
    }

    if (filters.area !== "all") {
      filtered = filtered.filter((row) => row.area.toLowerCase() === filters.area.toLowerCase())
    }

    setFilteredData(filtered)
  }

  // Chart 1: Prevalence of Overweight/Obese by State
  const obesityByState = filteredData
    .reduce((acc: any[], row) => {
      const existing = acc.find((item) => item.state === row.state)
      if (existing) {
        existing.women += row.women_overweight_obese_pct || 0
        existing.men += row.men_overweight_obese_pct || 0
        existing.count += 1
      } else {
        acc.push({
          state: row.state,
          women: row.women_overweight_obese_pct || 0,
          men: row.men_overweight_obese_pct || 0,
          count: 1,
        })
      }
      return acc
    }, [])
    .map((item) => ({
      state: item.state.length > 15 ? item.state.substring(0, 15) + "..." : item.state,
      women: (item.women / item.count).toFixed(1),
      men: (item.men / item.count).toFixed(1),
    }))
    .sort((a, b) => Number(b.women) - Number(a.women))
    .slice(0, 10)

  // Chart 2: Sample Distribution by Area
  const areaDistribution = filteredData.reduce((acc: any, row) => {
    const area = row.area || "Unknown"
    acc[area] = (acc[area] || 0) + 1
    return acc
  }, {})

  const pieData = Object.entries(areaDistribution).map(([area, count]) => ({
    name: area,
    value: count,
    percentage: (((count as number) / filteredData.length) * 100).toFixed(1),
  }))

  // Chart 3: Hypertension by Gender
  const hypertensionData = filteredData
    .reduce((acc: any[], row) => {
      const existing = acc.find((item) => item.state === row.state)
      if (existing) {
        existing.women += row.women_high_blood_pressure_pct || 0
        existing.men += row.men_high_blood_pressure_pct || 0
        existing.count += 1
      } else {
        acc.push({
          state: row.state,
          women: row.women_high_blood_pressure_pct || 0,
          men: row.men_high_blood_pressure_pct || 0,
          count: 1,
        })
      }
      return acc
    }, [])
    .map((item) => ({
      state: item.state.length > 12 ? item.state.substring(0, 12) + "..." : item.state,
      women: (item.women / item.count).toFixed(1),
      men: (item.men / item.count).toFixed(1),
    }))
    .sort((a, b) => Number(b.women) - Number(a.women))
    .slice(0, 8)

  // Chart 4: BMI Categories (simulated data based on overweight/obese percentages)
  const bmiData = [
    { category: "15-19 years", bmi: 21.2 },
    { category: "20-24 years", bmi: 22.8 },
    { category: "25-29 years", bmi: 24.1 },
    { category: "30-34 years", bmi: 25.3 },
    { category: "35-39 years", bmi: 25.8 },
    { category: "40-44 years", bmi: 26.2 },
    { category: "45-49 years", bmi: 26.0 },
  ]

  const COLORS = ["#06b6d4", "#ec4899", "#8b5cf6", "#10b981", "#f59e0b"]

  if (loading) {
    return (
      <div className="min-h-screen bg-slate-950 text-white p-6">
        <div className="container mx-auto">
          <div className="flex items-center justify-center h-64">
            <div className="text-cyan-400">Loading NFHS-5 data...</div>
          </div>
        </div>
      </div>
    )
  }

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="flex">
        {/* Sidebar */}
        <div className="fixed left-0 top-0 h-full p-6 bg-slate-950/50 backdrop-blur-sm border-r border-slate-700/50">
          <FilterSidebar onFiltersChange={handleFiltersChange} />
        </div>

        {/* Main Content */}
        <div className="flex-1 ml-96 p-6">
          <DashboardHeader
            title="NFHS-5 Data Explorer"
            subtitle="Exploratory data analysis of National Family Health Survey indicators"
          />

          {/* Visualization Grid */}
          <div className="grid lg:grid-cols-2 gap-8">
            {/* Chart 1: Obesity by State */}
            <GlassChart title="Prevalence of Overweight or Obese by State/UT">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={obesityByState} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                  <XAxis dataKey="state" stroke="#94a3b8" angle={-45} textAnchor="end" height={80} fontSize={10} />
                  <YAxis stroke="#94a3b8" />
                  <Tooltip
                    contentStyle={{
                      backgroundColor: "#1e293b",
                      border: "1px solid #06b6d4",
                      borderRadius: "8px",
                      color: "#fff",
                    }}
                  />
                  <Bar dataKey="women" fill="#06b6d4" name="Women" />
                  <Bar dataKey="men" fill="#ec4899" name="Men" />
                </BarChart>
              </ResponsiveContainer>
            </GlassChart>

            {/* Chart 2: Area Distribution */}
            <GlassChart title="Sample Distribution by Area">
              <ResponsiveContainer width="100%" height="100%">
                <PieChart>
                  <Pie
                    data={pieData}
                    cx="50%"
                    cy="50%"
                    innerRadius={60}
                    outerRadius={120}
                    paddingAngle={5}
                    dataKey="value"
                  >
                    {pieData.map((entry, index) => (
                      <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                    ))}
                  </Pie>
                  <Tooltip
                    contentStyle={{
                      backgroundColor: "#1e293b",
                      border: "1px solid #06b6d4",
                      borderRadius: "8px",
                      color: "#fff",
                    }}
                    formatter={(value: any, name: any, props: any) => [`${value} (${props.payload.percentage}%)`, name]}
                  />
                </PieChart>
              </ResponsiveContainer>
            </GlassChart>

            {/* Chart 3: Hypertension by Gender */}
            <GlassChart title="Prevalence of Hypertension by Gender">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={hypertensionData} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                  <XAxis dataKey="state" stroke="#94a3b8" angle={-45} textAnchor="end" height={80} fontSize={10} />
                  <YAxis stroke="#94a3b8" />
                  <Tooltip
                    contentStyle={{
                      backgroundColor: "#1e293b",
                      border: "1px solid #06b6d4",
                      borderRadius: "8px",
                      color: "#fff",
                    }}
                  />
                  <Bar dataKey="women" fill="#06b6d4" name="Women" />
                  <Bar dataKey="men" fill="#ec4899" name="Men" />
                </BarChart>
              </ResponsiveContainer>
            </GlassChart>

            {/* Chart 4: Mean BMI by Age Group */}
            <GlassChart title="Mean BMI by Age Group">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={bmiData} margin={{ top: 20, right: 30, left: 20, bottom: 60 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                  <XAxis dataKey="category" stroke="#94a3b8" angle={-45} textAnchor="end" height={80} fontSize={10} />
                  <YAxis stroke="#94a3b8" />
                  <Tooltip
                    contentStyle={{
                      backgroundColor: "#1e293b",
                      border: "1px solid #06b6d4",
                      borderRadius: "8px",
                      color: "#fff",
                    }}
                  />
                  <Bar dataKey="bmi" fill="#8b5cf6" />
                </BarChart>
              </ResponsiveContainer>
            </GlassChart>
          </div>
        </div>
      </div>
    </div>
  )
}
