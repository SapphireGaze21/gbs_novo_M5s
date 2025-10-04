import { GlassCard } from "./glass-card"
import { cn } from "@/lib/utils"

interface KpiCardProps {
  title: string
  value: string | number
  subtitle?: string
  className?: string
}

export function KpiCard({ title, value, subtitle, className }: KpiCardProps) {
  return (
    <GlassCard className={cn("p-6", className)}>
      <div className="text-center">
        <h3 className="text-sm font-medium text-slate-400 mb-2">{title}</h3>
        <div className="text-3xl font-bold text-white mb-1 bg-gradient-to-r from-cyan-400 to-teal-400 bg-clip-text text-transparent">
          {value}
        </div>
        {subtitle && <p className="text-xs text-slate-500">{subtitle}</p>}
      </div>
    </GlassCard>
  )
}
