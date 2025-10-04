import { GlassCard } from "./glass-card"
import { cn } from "@/lib/utils"
import type { ReactNode } from "react"

interface GlassChartProps {
  title: string
  children: ReactNode
  className?: string
}

export function GlassChart({ title, children, className }: GlassChartProps) {
  return (
    <GlassCard className={cn("p-6", className)}>
      <h3 className="text-lg font-semibold text-white mb-4">{title}</h3>
      <div className="w-full h-80">{children}</div>
    </GlassCard>
  )
}
