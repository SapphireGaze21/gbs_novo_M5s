import { cn } from "@/lib/utils"
import type { ReactNode } from "react"

interface GlassCardProps {
  children: ReactNode
  className?: string
}

export function GlassCard({ children, className }: GlassCardProps) {
  return (
    <div
      className={cn(
        "bg-slate-900/30 backdrop-blur-md border border-cyan-400/20 rounded-lg shadow-lg",
        "hover:border-cyan-400/40 transition-all duration-300",
        className,
      )}
    >
      {children}
    </div>
  )
}
