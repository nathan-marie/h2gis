package org.h2gis.functions.spatial.split;

import org.h2gis.api.DeterministicScalarFunction;
import org.locationtech.jts.geom.*;

import java.util.*;

public class ST_SubDivide extends DeterministicScalarFunction {

    static final GeometryFactory FACTORY = new GeometryFactory();

    public ST_SubDivide() {
        addProperty(PROP_REMARKS, "Divides geometry into parts using its internal envelope, " +
                "until each part can be represented using no more than max_vertices.\n" +
                "If no vertices apply a single recurse.");
    }

    @Override
    public String getJavaStaticMethod() {
        return "divide";
    }

    /**
     * Divide the geometry into quadrants (single pass, no vertex limit)
     */
    public static Geometry divide(Geometry geom) {
        return divideOnePass(geom);
    }

    /**
     * Divide the geometry recursively until each part has at most maxVertices.
     */
    public static Geometry divide(Geometry geom, int maxVertices) {
        if (geom == null || geom.isEmpty()) return geom;
        List<Geometry> parts = subdivideRecursive(geom, Math.max(0, maxVertices));
        if (parts == null || parts.isEmpty()) return geom.getFactory().createGeometryCollection();
        Geometry result = FACTORY.buildGeometry(parts);
        result.setSRID(geom.getSRID());
        return result;
    }

    // -------------------------------------------------------------------------
    // Core recursive subdivision
    // -------------------------------------------------------------------------

    /**
     * Recursively subdivides each sub-geometry of geom until no part exceeds maxVertices.
     */
    static List<Geometry> subdivideRecursive(Geometry geom, int maxVertices) {
        if (geom == null || geom.isEmpty()) return Collections.emptyList();

        Deque<Geometry> stack = new ArrayDeque<>();
        // Seed the stack with non-empty sub-geometries
        for (int i = 0; i < geom.getNumGeometries(); i++) {
            Geometry sub = geom.getGeometryN(i);
            if (!sub.isEmpty()) stack.push(sub);
        }

        List<Geometry> results = new ArrayList<>();
        while (!stack.isEmpty()) {
            Geometry slice = stack.pop();
            int nPts = vertexCount(slice);  // FIX: was using `geom` instead of `slice`
            if (nPts > maxVertices) {
                splitIntoQuadrants(slice, maxVertices, stack, results);
            } else {
                results.add(slice);
            }
        }
        return results;
    }

    // -------------------------------------------------------------------------
    // Single-pass subdivision (no vertex limit)
    // -------------------------------------------------------------------------

    private static Geometry divideOnePass(Geometry geom) {
        if (geom == null || geom.isEmpty()) return geom;

        List<Geometry> results = new ArrayList<>();
        for (int i = 0; i < geom.getNumGeometries(); i++) {
            Geometry sub = geom.getGeometryN(i);
            if (!sub.isEmpty() && (sub instanceof Polygon || sub instanceof LineString)) {
                quadrantIntersections(sub, results);
            } else if (!sub.isEmpty()) {
                results.add(sub);
            }
            // empty sub-geometries are silently dropped
        }

        Geometry res = FACTORY.buildGeometry(results);
        res.setSRID(geom.getSRID());
        return res;
    }

    // -------------------------------------------------------------------------
    // Shared geometry helpers
    // -------------------------------------------------------------------------

    /**
     * Returns the vertex count relevant for subdivision:
     *   - Polygon:    numPoints - 1  (closing point not counted)
     *   - LineString: numPoints
     *   - Other:      0              (never subdivided)
     */
    private static int vertexCount(Geometry g) {
        if (g instanceof Polygon)    return g.getNumPoints() - 1;
        if (g instanceof LineString) return g.getNumPoints();
        return 0;
    }

    /**
     * Splits a geometry into 2 or 4 quadrant intersections and routes each
     * result back to the stack (if too large) or to results (if small enough).
     */
    private static void splitIntoQuadrants(Geometry slice, int maxVertices,
                                           Deque<Geometry> stack, List<Geometry> results) {
        Envelope env = slice.getEnvelopeInternal();
        double minX = env.getMinX(), maxX = env.getMaxX(), midX = minX + (maxX - minX) / 2.0;
        double minY = env.getMinY(), maxY = env.getMaxY(), midY = minY + (maxY - minY) / 2.0;

        List<Envelope> quadrants;
        if (env.getHeight() == 0) {
            // Horizontal line: split left/right only
            quadrants = Arrays.asList(
                    new Envelope(minX, midX, midY, maxY),
                    new Envelope(midX, maxX, midY, maxY));
        } else if (env.getWidth() == 0) {
            // Vertical line: split top/bottom only
            quadrants = Arrays.asList(
                    new Envelope(minX, midX, minY, midY),
                    new Envelope(midX, maxX, minY, midY));
        } else {
            quadrants = Arrays.asList(
                    new Envelope(minX, midX, midY, maxY),
                    new Envelope(midX, maxX, midY, maxY),
                    new Envelope(minX, midX, minY, midY),
                    new Envelope(midX, maxX, minY, midY));
        }

        for (Envelope q : quadrants) {
            filterGeom(FACTORY.toGeometry(q).intersection(slice), maxVertices, stack, results);
        }
    }

    /**
     * Same quadrant logic for the single-pass variant (no size check).
     */
    private static void quadrantIntersections(Geometry slice, List<Geometry> results) {
        Envelope env = slice.getEnvelopeInternal();
        double minX = env.getMinX(), maxX = env.getMaxX(), midX = minX + (maxX - minX) / 2.0;
        double minY = env.getMinY(), maxY = env.getMaxY(), midY = minY + (maxY - minY) / 2.0;

        List<Envelope> quadrants;
        if (env.getHeight() == 0) {
            quadrants = Arrays.asList(
                    new Envelope(minX, midX, midY, maxY),
                    new Envelope(midX, maxX, midY, maxY));
        } else if (env.getWidth() == 0) {
            quadrants = Arrays.asList(
                    new Envelope(minX, midX, minY, midY),
                    new Envelope(midX, maxX, minY, midY));
        } else {
            quadrants = Arrays.asList(
                    new Envelope(minX, midX, midY, maxY),
                    new Envelope(midX, maxX, midY, maxY),
                    new Envelope(minX, midX, minY, midY),
                    new Envelope(midX, maxX, minY, midY));
        }

        for (Envelope q : quadrants) {
            Geometry inter = FACTORY.toGeometry(q).intersection(slice);
            if (!inter.isEmpty()) results.add(inter);
        }
    }

    /**
     * Dispatches each sub-geometry of a result back to the stack or results list.
     * Empty geometries and dimension-0 (point) results are discarded.
     */
    static void filterGeom(Geometry geom, int maxVertices,
                           Deque<Geometry> stack, List<Geometry> results) {
        for (int i = 0; i < geom.getNumGeometries(); i++) {
            Geometry sub = geom.getGeometryN(i);
            if (sub.isEmpty()) continue;  // FIX: drop empty sub-geometries
            int dim = sub.getDimension();
            if (dim == 2 || dim == 1) {
                int nPts = vertexCount(sub);
                if (nPts <= maxVertices) results.add(sub);
                else stack.push(sub);
            }
            // dim == 0 (points) silently dropped — cannot be subdivided
        }
    }
}