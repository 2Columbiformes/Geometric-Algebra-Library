using Godot;
using System;
using System.Linq;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Symbolics;
using Expr = MathNet.Symbolics.SymbolicExpression;


public static class LinearA {
	public static int[] rmElem(this int[] arr, int sel) {
		if (sel < 0 || sel >= arr.Length)
			throw new ArgumentOutOfRangeException(nameof(sel));
		int[] result = new int[arr.Length - 1];
		if (sel > 0)
			Array.Copy(arr, 0, result, 0, sel);
		if (sel < arr.Length - 1)
			Array.Copy(arr, sel + 1, result, sel, arr.Length - sel - 1);
		return result;
	}
	public static double[] Get(this double[,] matrix, int dim, int sel) {
		int n = matrix.GetLength(0);
		int m = matrix.GetLength(1);
		if (dim == 0) {
			if (sel < 0 || sel >= n) throw new ArgumentOutOfRangeException(nameof(sel));
			double[] res = new double[m];
			for (int j = 0; j < m; j++) res[j] = matrix[sel, j];
			return res;
		}
		else if (dim == 1) {
			if (sel < 0 || sel >= m) throw new ArgumentOutOfRangeException(nameof(sel));
			double[] res = new double[n];
			for (int i = 0; i < n; i++) res[i] = matrix[i, sel];
			return res;
		}
		else throw new ArgumentException("dim must be 0 (row) or 1 (column).");
	}
	public static int[] Get(this int[,] matrix, int dim, int sel) {
		int n = matrix.GetLength(0); 
		int m = matrix.GetLength(1);
		if (dim == 0) {
			if (sel < 0 || sel >= n) throw new ArgumentOutOfRangeException(nameof(sel));
			var res = new int[m];
			for (int j = 0; j < m; j++) res[j] = matrix[sel, j];
			return res;
		}
		else if (dim == 1) {
			if (sel < 0 || sel >= m) throw new ArgumentOutOfRangeException(nameof(sel));
			var res = new int[n];
			for (int i = 0; i < n; i++) res[i] = matrix[i, sel];
			return res;
		}
		else throw new ArgumentException("dim must be 0 (row) or 1 (column).");
	}
	public static double[,] Gets(this double[,] matrix, int dim, params int[] sel) {
		int n = matrix.GetLength(0);
		int m = matrix.GetLength(1);
		if (dim == 0) {
			n = sel.Length;
			foreach (int s in sel)
				if (s < 0 || s >= matrix.GetLength(0)) 
					throw new ArgumentOutOfRangeException(nameof(sel), "dim 0 index out of range."); 
			double[,] res = new double[n, m];
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					res[i, j] = matrix[sel[i], j];
			return res;
		}
		else if (dim == 1) {
			m = sel.Length;
			foreach (int s in sel)
				if (s < 0 || s >= matrix.GetLength(1)) 
					throw new ArgumentOutOfRangeException(nameof(sel), "dim 1 index out of range."); 
			double[,] res = new double[n, m];
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					res[i, j] = matrix[i, sel[j]];
			return res;
		}
		else throw new ArgumentException("dim must be 0 (row) or 1 (column).");
	}
	public static int[,] Gets(this int[,] matrix, int dim, params int[] sel) {
		int n = matrix.GetLength(0);
		int m = matrix.GetLength(1);
		if (dim == 0) {
			n = sel.Length;
			foreach (int s in sel)
				if (s < 0 || s >= matrix.GetLength(0)) 
					throw new ArgumentOutOfRangeException(nameof(sel), "dim 0 index out of range."); 
			int[,] res = new int[n, m];
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					res[i, j] = matrix[sel[i], j];
			return res;
		}
		else if (dim == 1) {
			m = sel.Length;
			foreach (int s in sel)
				if (s < 0 || s >= matrix.GetLength(1)) 
					throw new ArgumentOutOfRangeException(nameof(sel), "dim 1 index out of range."); 
			int[,] res = new int[n, m];
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					res[i, j] = matrix[i, sel[j]];
			return res;
		}
		else throw new ArgumentException("dim must be 0 (row) or 1 (column).");
	}
	public static double[] Sum(this double[,] v) {
		int n = v.GetLength(0); // rows
		int m = v.GetLength(1); // columns
		var result = new double[m];
		for (int j = 0; j < m; j++) {
			double sum = 0.0;
			for (int i = 0; i < n; i++)
				sum += v[i, j];
			result[j] = sum;
		}
		return result;
	}

	public static double Dot(this double[] v1, double[] v2) {
		if (v1.Length != v2.Length) throw new ArgumentException("Arrays must have the same length.");
		double result = 0.0;
		for (int i = 0; i < v1.Length; i++)
			result += v1[i] * v2[i];
		return result;
	}
	public static double Dot(this double[] v1) {
		double result = 0.0;
		for (int i = 0; i < v1.Length; i++)
			result += v1[i] * v1[i];
		return result;
	}
	public static double[] Norm(this double[] v1) {
		double mag = Math.Sqrt(v1.Dot());
		if (mag < 1e-12) throw new ArgumentException("Vector is too small to norm");
		for (int i = 0; i < v1.Length; i++)
			v1[i] = v1[i]/mag;
		return v1;
	}
	public static double[] Confine(this double[] v, double radius) {
		int n = v.Length;
		var res = new double[n];
		var len = Math.Sqrt(v.Dot());
		if (len > radius) {
			for (int i = 0; i < n; i++) {
				res[i] = v[i]/len*radius;
			}
			return res;
		}
		else return v;
	}
	public static double[,] Confine(this double[,] v, double radius) {
		int n = v.GetLength(0);
		int m = v.GetLength(1);
		var res = new double[n,m];
		for (int i = 0; i < n; i++) {
			double len = 0.0; 
			for (int j = 0; j < m; j++) 
				len += v[i, j] * v[i, j]; 
			len = Math.Sqrt(len);
			
			double scale = len > radius ? radius / len : 1.0;
			for (int j = 0; j < m; j++) res[i, j] = v[i, j] * scale;
		}
		return res;
	}
	public static double[] Cross(this double[] a, double[] b) {
		//if (a.Length != b.Length) throw new ArgumentException("Arrays must have the same length.");
		if (a.Length != 3) throw new ArgumentException("Arrays must have length 3.");
		return new double[] {
			a[1]*b[2] - a[2]*b[1],
			a[2]*b[0] - a[0]*b[2],
			a[0]*b[1] - a[1]*b[0]
		};
	}
	public static double[] Add(this double[] v1, double[] v2) {
		if (v1.Length != v2.Length) throw new ArgumentException("Arrays must have the same length.");
		double[] result = new double[v1.Length];
		for (int i = 0; i < v1.Length; i++)
			result[i] = v1[i] + v2[i];
		return result;
	}
	public static double[] Sub(this double[] v1, double[] v2) {
		if (v1.Length != v2.Length) throw new ArgumentException("Arrays must have the same length.");
		double[] result = new double[v1.Length];
		for (int i = 0; i < v1.Length; i++)
			result[i] = v1[i] - v2[i];
		return result;
	}
	public static double[] Mul(this double[] v, double factor) {
		double[] result = new double[v.Length];
		for (int i = 0; i < v.Length; i++)
			result[i] = v[i] * factor;
		return result;
	}
	public static double[] Div(this double[] v, double factor) {
		double[] result = new double[v.Length];
		for (int i = 0; i < v.Length; i++)
			result[i] = v[i] / factor;
		return result;
	}
	public static double Arg2pi(double x, double y) {
		// atan2 returns angle in [-π, π]
		double angle = Math.Atan2(y, x);
		// Convert to [0, 2π)
		if (angle < 0) angle += 2 * Math.PI;
		return angle;
	}

}
