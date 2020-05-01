"use strict";

class SpectralConformalParameterization {
	/**
	 * This class implements the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf spectral conformal parameterization} algorithm to flatten
	 * surface meshes with boundaries conformally.
	 * @constructor module:Projects.SpectralConformalParameterization
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the complex conformal energy matrix EC = ED - A.
	 * @private
	 * @method module:Projects.SpectralConformalParameterization#buildConformalEnergy
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	buildConformalEnergy() {
		
		let T = new ComplexTriplet(this.geometry.mesh.vertices.length, this.geometry.mesh.vertices.length);
		for (let bf of this.geometry.mesh.boundaries) {
			for (let h of bf.adjacentHalfedges()) {
				let vi = this.vertexIndex[h.vertex];
				let vj = this.vertexIndex[h.twin.vertex];
				T.addEntry(new Complex(0, 0.25), vi, vj);
				T.addEntry(new Complex(0, -0.25), vj, vi);
			}
		}
		let A = ComplexSparseMatrix.fromTriplet(T);
		return this.geometry.complexLaplaceMatrix(this.vertexIndex).timesComplex(new Complex(0.5)).minus(A);
	}

	/**
	 * Flattens the input surface mesh with 1 or more boundaries conformally.
	 * @method module:Projects.SpectralConformalParameterization#flatten
	 * @returns {Object} A dictionary mapping each vertex to a vector of planar coordinates.
	 */
	flatten() {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		let eigvec = Solvers.solveInversePowerMethod(this.buildConformalEnergy())
		let flattening = new Map();
		for (let v of vertices) {
			let c = eigvec.get(this.vertexIndex[v]);
			flattening[v] = new Vector(c.re, c.im);
		}
		
		// normalize flattening
		normalize(flattening, vertices);

		return flattening;
	}
}
