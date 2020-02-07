"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
		    let index = 0;
		    for (let v of mesh.vertices) {
			    v.index = index++;
		    }

		    index = 0;
		    for (let e of mesh.edges) {
			    e.index = index++;
		    }

		    index = 0;
		    for (let f of mesh.faces) {
			    f.index = index++;
		    }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
            let T = new Triplet(mesh.edges.length, mesh.vertices.length);
		    for (let e of mesh.edges) {
		        T.addEntry(1, e.index, e.halfedge.vertex.index)
		        T.addEntry(1, e.index, e.halfedge.twin.vertex.index)
		    }
            return SparseMatrix.fromTriplet(T);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
            let T = new Triplet(mesh.faces.length, mesh.edges.length);
            for (let f of mesh.faces) {
                T.addEntry(1, f.index, f.halfedge.edge.index)
                T.addEntry(1, f.index, f.halfedge.next.edge.index)
                T.addEntry(1, f.index, f.halfedge.next.next.edge.index)
            }
            return SparseMatrix.fromTriplet(T);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
            let u = DenseMatrix.zeros(this.mesh.vertices.length);
		    for (let v of subset.vertices) {
		        u.set(1, v);
		    }
		    return u;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
            let u = DenseMatrix.zeros(this.mesh.edges.length);
		    for (let e of subset.edges) {
		        u.set(1, e);
		    }
		    return u;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
            let u = DenseMatrix.zeros(this.mesh.faces.length);
            for (let f of subset.faces) {
                u.set(1, f);
            }
            return u;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
            let result = MeshSubset.deepCopy(subset);
            let es = this.A0.timesDense(this.buildVertexVector(result));
            for (let i = 0; i < es.nRows(); i++) {
                if (es.get(i) > 0) result.addEdge(i);
            }
            let fs = this.A1.timesDense(this.buildEdgeVector(result));
            for (let i = 0; i < fs.nRows(); i++) {
                if (fs.get(i) > 0) result.addFace(i);
            }
            return result;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
            let result = MeshSubset.deepCopy(subset);
            let es = this.A1.transpose().timesDense(this.buildFaceVector(result));
            for (let i = 0; i < es.nRows(); i++) {
                if (es.get(i) > 0) result.addEdge(i);
            }
            let vs = this.A0.transpose().timesDense(this.buildEdgeVector(result));
            for (let i = 0; i < vs.nRows(); i++) {
                if (vs.get(i) > 0) result.addVertex(i);
            }
            return result;
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
            let result = this.closure(this.star(subset));
            result.deleteSubset(this.star(this.closure(subset)));
            return result;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
          
            let vs = this.A0.transpose().timesDense(this.buildEdgeVector(subset));
            for (let i = 0; i < vs.nRows(); i++) {
                if (vs.get(i) > 0 && !subset.vertices.has(i)) {
                  return false;
                }
            }
            let es = this.A1.transpose().timesDense(this.buildFaceVector(subset));
            for (let i = 0; i < es.nRows(); i++) {
                if (es.get(i) > 0 && !subset.edges.has(i)) {
                  return false;
                }
            }
            return true;
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
            if (!this.isComplex(subset)) return -1;
            if (subset.faces.size > 0) {
                if (subset.equals(this.closure(new MeshSubset(new Set(), new Set(), subset.faces)))) return 2;
                else return -1;
            }
            else if (subset.edges.size > 0) {
                if (subset.equals(this.closure(new MeshSubset(new Set(), subset.edges)))) return 1;
                else return -1;
            }
            else return 0;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
            let d = this.isPureComplex(subset);
            let bd = new MeshSubset();
            if (d == 2) {
                let es = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                for (let i = 0; i < es.nRows(); i++) {
                    if (es.get(i) == 1) bd.addEdge(i);
                }
                return this.closure(bd);
            } else if (d == 1) {
                let vs = this.A0.transpose().timesDense(this.buildEdgeVector(subset));
                for (let i = 0; i < vs.nRows(); i++) {
                    if (vs.get(i) == 1) bd.addVertex(i);
                }
                return bd;
            } else if (d == 0) return bd;
        }
}
