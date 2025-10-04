"use client"

import { useState } from "react"
import { GlassCard } from "./glass-card"
import { Button } from "@/components/ui/button"
import { Label } from "@/components/ui/label"
import { Checkbox } from "@/components/ui/checkbox"
import { RadioGroup, RadioGroupItem } from "@/components/ui/radio-group"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { Slider } from "@/components/ui/slider"

interface FilterSidebarProps {
  onFiltersChange?: (filters: any) => void
}

export function FilterSidebar({ onFiltersChange }: FilterSidebarProps) {
  const [selectedStates, setSelectedStates] = useState<string[]>([])
  const [selectedArea, setSelectedArea] = useState<string>("all")
  const [ageRange, setAgeRange] = useState([15, 49])
  const [wealthQuintiles, setWealthQuintiles] = useState<string[]>([])

  const states = [
    "Andhra Pradesh",
    "Assam",
    "Bihar",
    "Gujarat",
    "Karnataka",
    "Kerala",
    "Maharashtra",
    "Tamil Nadu",
    "Uttar Pradesh",
    "West Bengal",
  ]

  const quintiles = ["Lowest", "Second", "Middle", "Fourth", "Highest"]

  const handleApplyFilters = () => {
    const filters = {
      states: selectedStates,
      area: selectedArea,
      ageRange,
      wealthQuintiles,
    }
    onFiltersChange?.(filters)
  }

  return (
    <div className="w-80 space-y-6">
      <GlassCard className="p-6">
        <h3 className="text-lg font-semibold text-white mb-4">Filters</h3>

        <div className="space-y-6">
          {/* State Selection */}
          <div>
            <Label className="text-slate-300 mb-3 block">State/UT</Label>
            <Select>
              <SelectTrigger className="bg-slate-800/50 border-slate-600 text-white">
                <SelectValue placeholder="Select states..." />
              </SelectTrigger>
              <SelectContent className="bg-slate-800 border-slate-600">
                {states.map((state) => (
                  <SelectItem key={state} value={state} className="text-white hover:bg-slate-700">
                    {state}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Area Selection */}
          <div>
            <Label className="text-slate-300 mb-3 block">Area</Label>
            <RadioGroup value={selectedArea} onValueChange={setSelectedArea}>
              <div className="flex items-center space-x-2">
                <RadioGroupItem value="all" id="all" className="border-slate-600 text-cyan-400" />
                <Label htmlFor="all" className="text-slate-300">
                  All
                </Label>
              </div>
              <div className="flex items-center space-x-2">
                <RadioGroupItem value="urban" id="urban" className="border-slate-600 text-cyan-400" />
                <Label htmlFor="urban" className="text-slate-300">
                  Urban
                </Label>
              </div>
              <div className="flex items-center space-x-2">
                <RadioGroupItem value="rural" id="rural" className="border-slate-600 text-cyan-400" />
                <Label htmlFor="rural" className="text-slate-300">
                  Rural
                </Label>
              </div>
            </RadioGroup>
          </div>

          {/* Age Range */}
          <div>
            <Label className="text-slate-300 mb-3 block">
              Age Range: {ageRange[0]} - {ageRange[1]} years
            </Label>
            <Slider value={ageRange} onValueChange={setAgeRange} max={65} min={15} step={5} className="w-full" />
          </div>

          {/* Wealth Quintiles */}
          <div>
            <Label className="text-slate-300 mb-3 block">Wealth Index Quintile</Label>
            <div className="space-y-2">
              {quintiles.map((quintile) => (
                <div key={quintile} className="flex items-center space-x-2">
                  <Checkbox
                    id={quintile}
                    checked={wealthQuintiles.includes(quintile)}
                    onCheckedChange={(checked) => {
                      if (checked) {
                        setWealthQuintiles([...wealthQuintiles, quintile])
                      } else {
                        setWealthQuintiles(wealthQuintiles.filter((q) => q !== quintile))
                      }
                    }}
                    className="border-slate-600 data-[state=checked]:bg-cyan-400 data-[state=checked]:border-cyan-400"
                  />
                  <Label htmlFor={quintile} className="text-slate-300">
                    {quintile}
                  </Label>
                </div>
              ))}
            </div>
          </div>

          <Button
            onClick={handleApplyFilters}
            className="w-full bg-gradient-to-r from-cyan-500 to-teal-500 hover:from-cyan-600 hover:to-teal-600 text-white"
          >
            Apply Filters
          </Button>
        </div>
      </GlassCard>
    </div>
  )
}
