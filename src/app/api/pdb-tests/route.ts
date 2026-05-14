import { db } from "@/lib/db";
import { NextRequest, NextResponse } from "next/server";

// GET /api/pdb-tests — list all tests
export async function GET() {
  try {
    const tests = await db.pdbTest.findMany({
      orderBy: { updatedAt: "desc" },
    });
    return NextResponse.json(tests);
  } catch (error) {
    console.error("Failed to list PDB tests:", error);
    return NextResponse.json({ error: "Failed to list tests" }, { status: 500 });
  }
}

// POST /api/pdb-tests — create a new test
export async function POST(request: NextRequest) {
  try {
    const body = await request.json();
    const { pdbId: inputPdbId } = body;

    // Generate a unique test ID like PDB-001, PDB-002, ...
    // Query the max numeric suffix from existing testIds to avoid conflicts after deletions
    let nextNum = 1;
    const allTests = await db.pdbTest.findMany({ select: { testId: true } });
    for (const t of allTests) {
      const match = t.testId.match(/^PDB-(\d+)$/);
      if (match) {
        const num = parseInt(match[1], 10);
        if (num >= nextNum) nextNum = num + 1;
      }
    }
    const testId = `PDB-${String(nextNum).padStart(3, "0")}`;

    const test = await db.pdbTest.create({
      data: {
        testId,
        pdbId: inputPdbId || null,
        state: "edit",
        currentStep: 1,
      },
    });
    return NextResponse.json(test);
  } catch (error) {
    console.error("Failed to create PDB test:", error);
    return NextResponse.json({ error: "Failed to create test" }, { status: 500 });
  }
}
