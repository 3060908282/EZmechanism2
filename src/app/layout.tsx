import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";
import { Toaster } from "@/components/ui/toaster";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: "M-CSA Mechanism Predictor | Enzyme Reaction Mechanism Prediction",
  description:
    "Predict enzyme reaction mechanisms using M-CSA curated catalytic mechanism rules. Input substrate SMILES and discover potential reaction pathways.",
  keywords: [
    "M-CSA",
    "enzyme",
    "mechanism",
    "prediction",
    "reaction",
    "SMARTS",
    "catalysis",
    "biochemistry",
  ],
  icons: {
    icon: "/molecule-icon.svg",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased bg-background text-foreground`}
      >
        {/* 3Dmol.js — loaded as external UMD script to avoid large Turbopack chunks */}
        {/* Pre-initialize window["3Dmol"] so internal modules (e.g. workerString setup) don't crash */}
        <script
          dangerouslySetInnerHTML={{
            __html: `window["3Dmol"] = window["3Dmol"] || {};`,
          }}
        />
        <script
          src="/3Dmol.min.js"
          defer
        />

        {children}
        <Toaster />
      </body>
    </html>
  );
}

