import dynamic from "next/dynamic"

const MedicalNLPAnalysis = dynamic(
  () => import("@/components/medical-nlp-analysis").then((m: any) => m.default ?? m.MedicalNLPAnalysis),
  { ssr: false },
)

export default function NLPInsightsPage() {
  return (
    <main className="px-4 py-6 md:px-8">
      <MedicalNLPAnalysis />
    </main>
  )
}
