"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		
    let Vs = geometry.mesh.vertices;
    let T = new Triplet(Vs.length, Vs.length);
    for (let v of Vs) {
      T.addEntry(geometry.barycentricDualArea(v), vertexIndex[v], vertexIndex[v]);
    }
		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
    let Es = geometry.mesh.edges;
    let T = new Triplet(Es.length, Es.length);
    for (let e of Es) {
      T.addEntry(0.5 * (geometry.cotan(e.halfedge) + geometry.cotan(e.halfedge.twin)), edgeIndex[e], edgeIndex[e]);
    }
		return SparseMatrix.fromTriplet(T);
    
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
    let Fs = geometry.mesh.faces;
    let T = new Triplet(Fs.length, Fs.length);
    for (let f of Fs) {
      T.addEntry(1.0 / geometry.area(f), faceIndex[f], faceIndex[f]);
    }
		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		
    let Vs = geometry.mesh.vertices;
    let Es = geometry.mesh.edges;
    
    let T = new Triplet(Es.length, Vs.length);
    for (let e of Es) {
      T.addEntry(-1, edgeIndex[e], vertexIndex[e.halfedge.vertex]);
      T.addEntry(1, edgeIndex[e], vertexIndex[e.halfedge.twin.vertex]);
    }
    
		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
						
    let Es = geometry.mesh.edges;
    let Fs = geometry.mesh.faces;
    
    let T = new Triplet(Fs.length, Es.length);
    for (let f of Fs) {
      let h = f.halfedge;
      for (let i = 0; i < 3; i++) {
        let e = h.edge;
        let dir = 1;
        if (e.halfedge != h) {
          dir = -1;
        }
        T.addEntry(dir, faceIndex[f], edgeIndex[e]);
        h = h.next;
      }
    }
		return SparseMatrix.fromTriplet(T);
    
	}
}
