"use client";

import { Card, CardContent, CardHeader } from "@/components/ui/card";
import { Skeleton } from "@/components/ui/skeleton";

export default function LoadingSkeleton() {
  return (
    <Card className="border-emerald-100 shadow-sm">
      <CardHeader>
        <Skeleton className="h-5 w-28" />
        <Skeleton className="h-3 w-48 mt-1.5" />
      </CardHeader>
      <CardContent className="space-y-3">
        <Skeleton className="h-8 w-full rounded-lg" />
        <div className="space-y-2">
          <Skeleton className="h-16 w-full" />
          <Skeleton className="h-16 w-full" />
          <Skeleton className="h-16 w-full" />
        </div>
        <div className="grid grid-cols-3 gap-3 mt-3">
          <Skeleton className="h-14 rounded-lg" />
          <Skeleton className="h-14 rounded-lg" />
          <Skeleton className="h-14 rounded-lg" />
        </div>
      </CardContent>
    </Card>
  );
}
