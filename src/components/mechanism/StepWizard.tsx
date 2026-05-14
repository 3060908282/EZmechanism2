"use client";

import React from "react";
import { CheckCircle2 } from "lucide-react";
import { STEPS } from "@/lib/constants";

export default function StepWizard({ currentStep, onStepClick }: { currentStep: number; onStepClick: (step: number) => void }) {
  return (
    <div className="w-full py-4 px-2">
      <div className="flex items-center justify-between max-w-2xl mx-auto">
        {STEPS.map((step, idx) => {
          const isCompleted = step.id < currentStep;
          const isActive = step.id === currentStep;
          const isFuture = step.id > currentStep;
          const canClick = !isFuture;
          return (
            <React.Fragment key={step.id}>
              <button
                onClick={() => canClick && onStepClick(step.id)}
                disabled={isFuture}
                className={`flex flex-col items-center gap-1 group transition-all ${
                  isFuture ? "opacity-50 cursor-not-allowed" : "cursor-pointer"
                }`}
              >
                <div
                  className={`w-9 h-9 rounded-full flex items-center justify-center text-sm font-bold transition-all border-2 ${
                    isCompleted
                      ? "bg-emerald-500 border-emerald-500 text-white shadow-md group-hover:shadow-lg"
                      : isActive
                        ? "bg-violet-600 border-violet-600 text-white shadow-lg ring-4 ring-violet-100"
                        : "bg-white border-gray-200 text-gray-400 group-hover:border-gray-300"
                  }`}
                >
                  {isCompleted ? <CheckCircle2 className="w-4 h-4" /> : <span>{step.id}</span>}
                </div>
                <span
                  className={`text-xs font-medium text-center max-w-[70px] leading-tight ${
                    isActive ? "text-violet-700" : isCompleted ? "text-emerald-600" : "text-gray-400"
                  }`}
                >
                  {step.label}
                </span>
              </button>
              {idx < STEPS.length - 1 && (
                <div
                  className={`flex-1 h-0.5 mx-1 sm:mx-2 rounded-full transition-all ${
                    step.id < currentStep ? "bg-emerald-400" : "bg-gray-200"
                  }`}
                />
              )}
            </React.Fragment>
          );
        })}
      </div>
    </div>
  );
}
