import { db } from "@/lib/db";
import { NextRequest, NextResponse } from "next/server";

// GET /api/pdb-tests/[id] — get single test
export async function GET(
  _request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id } = await params;
    const test = await db.pdbTest.findUnique({ where: { id } });
    if (!test) {
      return NextResponse.json({ error: "Test not found" }, { status: 404 });
    }
    return NextResponse.json(test);
  } catch (error) {
    console.error("Failed to get PDB test:", error);
    return NextResponse.json({ error: "Failed to get test" }, { status: 500 });
  }
}

// PUT /api/pdb-tests/[id] — update test (auto-save)
export async function PUT(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id } = await params;
    const body = await request.json();

    // Remove fields that should not be updated from client
    const { id: _id, testId: _testId, createdAt: _createdAt, updatedAt: _updatedAt, ...data } = body as Record<string, unknown>;

    const test = await db.pdbTest.update({
      where: { id },
      data,
    });
    return NextResponse.json(test);
  } catch (error) {
    console.error("Failed to update PDB test:", error);
    return NextResponse.json({ error: "Failed to update test" }, { status: 500 });
  }
}

// DELETE /api/pdb-tests/[id] — delete test
export async function DELETE(
  _request: NextRequest,
  { params }: { params: Promise<{ id: string }> }
) {
  try {
    const { id } = await params;
    await db.pdbTest.delete({ where: { id } });
    return NextResponse.json({ success: true });
  } catch (error) {
    console.error("Failed to delete PDB test:", error);
    return NextResponse.json({ error: "Failed to delete test" }, { status: 500 });
  }
}
