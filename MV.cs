using Godot;
using System;
using System.Linq;
using System.Collections.Generic;
using MathNet.Symbolics;
using System.Numerics;
using Expr = MathNet.Symbolics.SymbolicExpression;


// Multivector in Cl(p,q)
public static class GA {
	public static int n;                // dimension
	public static int size;             // 2^n
	public static int[] signature;      // length n, entries +1 or -1
	public static int[] flatGrade;
	public static int[][] gradeMap;     // grade -> indices
	public static int[,] mask;          // product blade index
	public static int[,] sign;          // product sign
	public static int[] bladeSign;		 //convention sign of multivector basis, designated by hodge dual
	// --- Initialization ---
	static int[] BuildFlatGradeMap(int n) {
		int size = 1 << n; // 2^n total blades
		var ordered = new List<int>(size);
		// For each grade in order
		for (int g = 0; g <= n; g++) {
			for (int i = 0; i < size; i++) {
				if (BitOperations.PopCount((uint)i) == g)
					ordered.Add(i);
			}
		}
		return ordered.ToArray();
	}
	static int[][] BuildGradeMapFromFlat(int n, int[] flat) {
		// Infer n from flat length (since flat.Length = 2^n)
		var result = new int[n + 1][];
		int index = 0;
		for (int g = 0; g <= n; g++) {
			// number of blades of grade g = C(n, g)
			int count = Choose(n, g);
			result[g] = new int[count];
			Array.Copy(flat, index, result[g], 0, count);
			index += count;
		}
		return result;
	}
	// Inline binomial coefficient (no separate helper needed elsewhere)
	static int Choose(int n, int k) {
		int res = 1;
		for (int i = 1; i <= k; i++) 
			res = res * (n - (k - i)) / i;
		return res;
	}
	public static void Init(int p, int q, bool swap = false) {
		n = p + q;
		size = 1 << n;
		if (swap)
			signature = Enumerable.Repeat(-1, p).Concat(Enumerable.Repeat(+1, q)).ToArray();
		else
			signature = Enumerable.Repeat(+1, p).Concat(Enumerable.Repeat(-1, q)).ToArray();
		// build raw grade map and flatten
		flatGrade = BuildFlatGradeMap(n);
		gradeMap = BuildGradeMapFromFlat(n, flatGrade);
		//flat = gradeMap.SelectMany(g => g).ToArray(); // keep flat as a field if needed
		mask = new int[size, size];
		sign = new int[size, size];
		bladeSign = Enumerable.Repeat(1, size).ToArray(); // initialize all to +1
		// build multiplication tables
		for (int I = 0; I < size; I++)
			for (int J = 0; J < size; J++) {
				mask[I, J] = I ^ J; //exclusion of bit representations
				sign[I, J] = SwapSign(I, J) * MetricSign(I, J, signature);
			}
	}
	public static void ReorderGrade(int grade, params int[] order) {
		// Restrictions: grade 0 and grades > n/2 cannot be reordered
		if (grade == 0 || grade > n / 2)
			throw new InvalidOperationException($"Grade {grade} cannot be reordered.");
		var indices = gradeMap[grade];
		int len = indices.Length;
		if (order == null || order.Length != len)
			throw new ArgumentException("array length must match # of blades in grade.");
		// Validate that order contains all values from 0 to len - 1 exactly once
		var seen = new bool[len];
		foreach (int pos in order) {
			if (pos < 0 || pos >= len)
				throw new ArgumentException($"Order array contains out-of-bounds index: {pos}");
			if (seen[pos])
				throw new ArgumentException($"Order array contains duplicate index: {pos}");
			seen[pos] = true;
		}
		// Apply permutation
		var reordered = new int[len];
		for (int i = 0; i < len; i++)
			reordered[i] = indices[order[i]];
		gradeMap[grade] = reordered;
		var flatList = new List<int>();
		for (int g = 0; g <= n; g++)
			flatList.AddRange(gradeMap[g]);
		flatGrade = flatList.ToArray();
	}
	public static void BladeSignHodge() {
		int full = (1 << n) - 1;
		var flat = flatGrade;
		for (int i = 0; i < size / 2; i++) {
			int I = flat[i], D = full ^ I;
			bladeSign[I] = 1;
			bladeSign[D] = SwapSign(I, D);
			flat[i + size / 2] = D;
		}
		// rebuild gradeMap from reordered flat
		gradeMap = Enumerable.Range(0, n + 1)
			.Select(g => flat.Where(idx => BitOperations.PopCount((uint)idx) == g).ToArray())
			.ToArray();
		flatGrade = gradeMap.SelectMany(g => g).ToArray();
		for (int I = 0; I < size; I++) {
			for (int J = 0; J < size; J++) {
				int K = I ^ J;
				// base sign from swaps + metric
				int baseSign = SwapSign(I, J) * MetricSign(I, J, signature);
				// adjust by bladeSigns
				sign[I, J] = baseSign * bladeSign[I] * bladeSign[J] * bladeSign[K];
			}
		}
	}
	// --- Sign helpers ---
	public static void FlipBladeSign(int bladeMask) {
		bladeSign[bladeMask] *= -1;
		// Rebuild multiplication tables with new convention
		for (int I = 0; I < size; I++)
			for (int J = 0; J < size; J++) {
				int K = I ^ J;
				sign[I, J] *= bladeSign[I] * bladeSign[J] * bladeSign[K];
			}
	}
	public static int SwapSign(int I, int J) {
		int swaps = 0;
		for (int b = 0; b < n; b++) {
			if (((J >> b) & 1) == 0) continue;
			for (int a = b + 1; a < n; a++)
				swaps += (I >> a) & 1;
		}
		return (swaps & 1) == 0 ? +1 : -1;
	}
	private static int MetricSign(int I, int J, int[] signature) {
		int s = 1;
		int common = I & J; //common bits
		for (int b = 0; b < signature.Length; b++)
			if (((common >> b) & 1) != 0) //if bit b is 1 in int common
				s *= signature[b];
		return s;
	}
	// --- Pretty print ---
	public static void PrintDebugInfo() {
		GD.Print($"n = {GA.n}");
		GD.Print($"GA.size = {GA.size}");
		GD.Print("signature = [" + string.Join(", ", GA.signature) + "]");
		GD.Print("flatGrade = [" + string.Join(", ", GA.flatGrade) + "]");
		GD.Print("gradeMap:");
		for (int g = 0; g < GA.gradeMap.Length; g++) {
			GD.Print($"  grade {g}: [" + string.Join(", ", GA.gradeMap[g]) + "]");
		}
		GD.Print("GA.mask table:");
		for (int i = 0; i < GA.size; i++) {
			string row = "";
			for (int j = 0; j < GA.size; j++) {
				row += GA.mask[i, j].ToString().PadLeft(3);
				if (j < GA.size - 1) row += " ";
			}
			GD.Print(row);
		}

		GD.Print("sign table:");
		for (int i = 0; i < GA.size; i++) {
			string row = "";
			for (int j = 0; j < GA.size; j++) {
				row += GA.sign[i, j].ToString().PadLeft(3);
				if (j < GA.size - 1) row += " ";
			}
			GD.Print(row);
		}

		GD.Print("GA.bladeSign = [" + string.Join(", ", GA.bladeSign) + "]");
	}
}

//numeric multivectors
public struct MVd {
	public double[] c;   // coefficients, length = GA.size
	// --- Constructors ---
	public MVd(params double[] components) {
		if (components.Length != GA.size)
			throw new ArgumentException($"MVd requires exactly {GA.size} components for Cl({GA.n}).");
		c = new double[GA.size];
		Array.Copy(components, c, GA.size);
	}
	public MVd(int grade, params double[] comps) {
		var indices = GA.gradeMap[grade];
		if (comps.Length != indices.Length)
			throw new ArgumentException($"Grade {grade} requires {indices.Length} component(s).");
		c = new double[GA.size];
		for (int i = 0; i < indices.Length; i++)
			c[indices[i]] = comps[i];
	}
	public static MVd CreateLie(int grade, params double[] comps) {
		if (grade < 1) throw new ArgumentException($"Grade must be larger than 0");
		var indices = GA.gradeMap[grade];
		if (comps.Length != indices.Length + 1)
			throw new ArgumentException($"Grade {grade} + scaler requires {indices.Length+1} component(s).");
		var res = new double[GA.size];
		for (int i = 0; i < indices.Length; i++){
			res[indices[i]] = comps[i+1];
		}
		res[indices[0]] = comps[0];
		return new MVd(res);
	}
	public static MVd Zero => new MVd(new double[GA.size]);

	public static MVd basis(int blade) {
		var mv = Zero;
		mv.c[blade] = 1.0;
		return mv;
	}

	// --- Geometric Product ---
	public static MVd operator *(MVd A, MVd B) {
		var r = new double[GA.size];
		for (int I = 0; I < GA.size; I++) {
			double a = A.c[I];
			if (a == 0) continue;
			for (int J = 0; J < GA.size; J++) {
				double b = B.c[J];
				if (b == 0) continue;
				int K = GA.mask[I, J];
				int sgn = GA.sign[I, J];
				r[K] += sgn * a * b;
			}
		}
		return new MVd(r);
	}

	// --- Dot product (scalar inner product) ---
	public double Dot(MVd B) {
		double sum = 0;
		for (int i = 0; i < GA.size; i++)
			sum += this.c[i] * B.c[i];
		return sum;
	}

	public MVd HodgeDual() {
		var r = new double[GA.size];
		int full = (1 << GA.n) - 1;
		for (int I = 0; I < GA.size; I++) {
			double a = c[I];
			if (a == 0) continue;
			int D = full ^ I;
			int sgn = GA.SwapSign(I, D);
			sgn *= GA.bladeSign[I] * GA.bladeSign[D];
			r[D] += sgn * a;
		}
		return new MVd(r);
	}

	// --- Left Contraction ---
	public MVd LInner(MVd B) {
		var r = new double[GA.size];
		for (int I = 0; I < GA.size; I++) {
			double a = c[I];
			if (a == 0) continue;
			int gI = BitOperations.PopCount((uint)I);
			for (int J = 0; J < GA.size; J++) {
				double b = B.c[J];
				if (b == 0) continue;
				int gJ = BitOperations.PopCount((uint)J);
				if (gI > gJ) continue;
				int K = GA.mask[I, J];
				int gK = BitOperations.PopCount((uint)K);
				if (gK != gJ - gI) continue;
				int sgn = GA.sign[I, J];
				r[K] += sgn * a * b;
			}
		}
		return new MVd(r);
	}
	public MVd Outer(MVd B) {
		var r = new double[GA.size];
		for (int I = 0; I < GA.size; I++) {
			double a = c[I];
			if (a == 0) continue;
			int gI = BitOperations.PopCount((uint)I);
			for (int J = 0; J < GA.size; J++) {
				double b = B.c[J];
				if (b == 0) continue;
				int gJ = BitOperations.PopCount((uint)J);
				int K = GA.mask[I, J];
				int gK = BitOperations.PopCount((uint)K);
				if (gK == gI + gJ) {
					int sgn = GA.sign[I, J];
					r[K] += sgn * a * b;
				}
			}
		}
		return new MVd(r);
	}
	
	public double[] Grade(int grade) {
		var idx = GA.gradeMap[grade];
		var comps = new double[idx.Length];
		for (int i = 0; i < idx.Length; i++)
			comps[i] = c[idx[i]];
		return comps;
	}
	
	public MVd Conj() {
		var r = new double[GA.size];
		for (int I = 0; I < GA.size; I++) {
			int g = BitOperations.PopCount((uint)I);
			int sgn = ((g * (g - 1) / 2) % 2 == 0) ? +1 : -1;
			r[I] = sgn * c[I];
		}
		return new MVd(r);
	}
	public bool IsInvertible() {
		var rev = this.Conj();
		var prod = this * rev;
		// Check that all non-scalar components are zero
		for (int i = 1; i < GA.size; i++) {
			if (prod.c[i] != 0.0)
				return false;   // not purely scalar → no simple inverse
		}
		// Scalar part must be non-zero
		return prod.c[0] != 0.0;
	}
	// --- Inverse ---
	public MVd Inv() {
		var rev = this.Conj();
		var prod = this * rev;
		// Check scalar-only
		for (int i = 1; i < GA.size; i++) {
			if (prod.c[i] != 0.0)
				throw new InvalidOperationException("Multivector is not a versor; inverse not defined by M~ / (M M~).");
		}
		double norm = prod.c[0];
		if (norm == 0.0)
			throw new InvalidOperationException("Multivector not invertible (zero norm).");
		return (1.0 / norm) * rev;
	}

	public static MVd Ave(params MVd[] mvs){
		MVd res = MVd.Zero;
		int n = mvs.Length;
		for (int i = 0; i < n; i++)
			res += mvs[i];
		return res / n;
	}
	public static MVd OuterProd(double[,] points, int[] selection) {
		if (selection == null || selection.Length < 2) return MVd.Zero;
		int n = points.GetLength(1);
		if (n != GA.n) throw new ArgumentException($"Point dimension {n} does not match GA dimension {GA.n}.");
		// First edge: p[sel[1]] - p[sel[0]]
		MVd result = new MVd(1, points.Get(0, selection[0]));
		// Wedge chain
		for (int k = 1; k < selection.Length; k++) {
			MVd edge = new MVd(1, points.Get(0, selection[k]));
			result = result.Outer(edge);
		}
		return result;
	}
	public static MVd Measure(double[,] points, int[] selection) {
		if (selection == null || selection.Length < 2) return MVd.Zero;
		int n = points.GetLength(1);
		if (n != GA.n) throw new ArgumentException($"Point dimension {n} does not match GA dimension {GA.n}.");
		// First edge: p[sel[1]] - p[sel[0]]
		MVd result = new MVd(1, points.Get(0, selection[1]))
				   - new MVd(1, points.Get(0, selection[0]));
		// Wedge chain
		for (int k = 2; k < selection.Length; k++) {
			MVd edge = new MVd(1, points.Get(0, selection[k]))
					 - new MVd(1, points.Get(0, selection[k - 1]));
			result = result.Outer(edge);
		}
		return result;
	}
	public static double[] Circ(double[,] points, int[] sel) {
		int m = sel.Length; // number of points
		int n = points.GetLength(1); // dimension of each point
		if (n != GA.n) throw new ArgumentException($"Point dimension {n} does not match GA dimension {GA.n}.");
		if (m <= 1) return new double[n];
		MVd result = MVd.Zero;
		if (m >= 3){
			for (int i = 0; i < m; i++) {
				result -= ((1.0-(2.0*(i%2.0)))/2.0) * points.Get(0,sel[i]).Dot() * Measure(points, sel.rmElem(i));
			}
			if (m <= n) {
				result += OuterProd(points, sel); //same as Measure_n+1(points, 0)
			}
		}
		else {
			var p1 = points.Get(0,sel[0]);
			var p2 = points.Get(0,sel[1]);
			var p1m = new MVd(1,p1);
			var p2m = new MVd(1,p2);
			result += p1.Dot()*p2m - p1m*p1.Dot(p2) - p2.Dot()*p1m + p2m*p2.Dot(p1); 
			result = result/(p1m.Outer(p2m)*2);
			return result.Grade(1);
		}
		result = result/Measure(points, sel); //GA invertsion
		return result.Grade(1);
	}
	
	// --- Operators ---
	public static MVd operator +(MVd a, MVd b) {
		var r = new double[GA.size];
		for (int i = 0; i < GA.size; i++)
			r[i] = a.c[i] + b.c[i];
		return new MVd(r);
	}
	public static MVd operator -(MVd a, MVd b) {
		var r = new double[GA.size];
		for (int i = 0; i < GA.size; i++)
			r[i] = a.c[i] - b.c[i];
		return new MVd(r);
	}
	public static MVd operator -(MVd a) {
		return (-1)*a;
	}
	public static MVd operator *(double s, MVd a) {
		var r = new double[GA.size];
		for (int i = 0; i < GA.size; i++)
			r[i] = s * a.c[i];
		return new MVd(r);
	}

	public static MVd operator *(MVd a, double s) => s * a;
	public static MVd operator /(MVd a, double s) {
		var r = new double[GA.size];
		double inv = 1.0 / s;
		for (int i = 0; i < GA.size; i++)
			r[i] = a.c[i] * inv;
		return new MVd(r);
	}
	public static MVd operator /(MVd A, MVd B) {
		return A * B.Inv();
	}
	public override string ToString() {
		return this.ToString(1);
	}
	public string ToString(int startIndex) {
		var sb = new System.Text.StringBuilder();
		bool first = true;
		// Flatten GA.gradeMap into dual-aware order
		var flat = GA.gradeMap.SelectMany(g => g).ToArray();
		foreach (var idx in flat) {
			double simple = c[idx];
			if (simple == 0.0) continue;
			if (!first)
				sb.Append(" + \n");
			sb.Append("[");
			sb.Append(simple.ToString("G"));   // numeric formatting
			sb.Append("]");
			if (idx != 0) {
				// print explicit blade indices
				string mask = Convert.ToString(idx, 2).PadLeft(GA.n, '0');
				char[] chars = mask.ToCharArray();
				Array.Reverse(chars);
				var blade = new System.Text.StringBuilder();
				for (int bit = 0; bit < chars.Length; bit++) {
					if (chars[bit] == '1')
						blade.Append(bit + startIndex);
				}
				if (blade.Length > 0) {
					if (GA.bladeSign[idx] < 0) sb.Append("*-e");
					else sb.Append("*e");
					sb.Append(blade.ToString());
				}
			}
			first = false;
		}
		return first ? "0" : sb.ToString();
	}
	
}

//Symbolic multivectors 
public struct MV {
	public Expr[] c; // coefficients, length = 2^n
	// static context for this algebra
	// --- Constructors ---
	public MV(params Expr[] components) {
		c = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++)
			c[i] = Expr.Zero;
		if (components.Length != GA.size)
			throw new ArgumentException($"MV requires exactly {GA.size} components for Cl({GA.n}).");
		for (int i = 0; i < GA.size; i++)
			c[i] = components[i] ?? Expr.Zero;
	}
	public static MV Zero => new MV();
	public MV(string prefix, int grade, int start = 1) {
		c = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++)
			c[i] = Expr.Zero;
		for (int i = 0; i < GA.gradeMap[grade].Length; i++) {
			c[GA.gradeMap[grade][i]] = Expr.Variable($"{prefix}{(start + i).ToString()}");
		}
	}
	public MV(int grade, params Expr[] comps) {
		c = new Expr[GA.size];
		var indices = GA.gradeMap[grade]; // already dual-aware ordering
		if (comps.Length != indices.Length)
			throw new ArgumentException($"Grade {grade} requires {indices.Length} component(s).");
		var arr = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) arr[i] = Expr.Zero;
		for (int i = 0; i < indices.Length; i++)
			arr[indices[i]] = comps[i];
		c = arr;
	}
	public static MV CreateLie(string prefix, int grade, int start = 0) {
		var res = new Expr[GA.size];
		var indices = GA.gradeMap[grade];
		for (int i = 0; i < indices.Length; i++) {
			res[indices[i]] = Expr.Variable($"{prefix}{(start + i + 1).ToString()}");
		}
		return new MV(0, Expr.Variable($"{prefix}{start.ToString()}")) + new MV(res);
	}
	
	//GA functions
	public static MV operator *(MV A, MV B) {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) r[i] = Expr.Zero;
		for (int I = 0; I < GA.size; I++) {
			if (A.c[I].IsZero()) continue;
			for (int J = 0; J < GA.size; J++) {
				if (B.c[J].IsZero()) continue;
				int K = GA.mask[I, J];
				int sgn = GA.sign[I, J];
				r[K] += Expr.Integer(sgn) * A.c[I] * B.c[J];
			}
		}
		return new MV(r);
	}

	public Expr Dot(MV b) {
		Expr result = Expr.Zero;
		for (int i = 0; i < GA.size; i++)
			result += this.c[i] * b.c[i];
		return result;
	}

	
	public MV HodgeDual() {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) r[i] = Expr.Zero;
		int full = (1 << GA.n) - 1; // GA.mask of e0e1...en-1
		for (int I = 0; I < GA.size; I++) {
			var a = c[I];
			if (a.IsZero()) continue;
			int D = full ^ I;         // complement GA.mask
			int sgn = GA.SwapSign(I, D); // permutation GA.sign
			// Apply GA.bladeSign convention
			sgn *= GA.bladeSign[I] * GA.bladeSign[D];
			r[D] += Expr.Integer(sgn) * a;
		}
		return new MV(r);
	}
	
	public MV LInner(MV B) {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) r[i] = Expr.Zero;
		for (int I = 0; I < GA.size; I++) {
			if (c[I].IsZero()) continue;
			int gI = System.Numerics.BitOperations.PopCount((uint)I);
			for (int J = 0; J < GA.size; J++) {
				if (B.c[J].IsZero()) continue;
				int gJ = System.Numerics.BitOperations.PopCount((uint)J);
				if (gI > gJ) continue; // left contraction requires r <= s
				int K = GA.mask[I, J];
				int gK = System.Numerics.BitOperations.PopCount((uint)K);
				if (gK != gJ - gI) continue; // project onto grade s - r
				int sgn = GA.sign[I, J];
				r[K] += Expr.Integer(sgn) * c[I] * B.c[J];
			}
		}
		return new MV(r);
	}
	public MV Outer(MV B) {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) r[i] = Expr.Zero;
		for (int I = 0; I < GA.size; I++) {
			if (c[I].IsZero()) continue;
			int gI = BitOperations.PopCount((uint)I);
			for (int J = 0; J < GA.size; J++) {
				if (B.c[J].IsZero()) continue;
				int gJ = BitOperations.PopCount((uint)J);
				int K = GA.mask[I, J];
				int gK = BitOperations.PopCount((uint)K);
				// Outer product only keeps terms where grades add
				if (gK == gI + gJ) {
					int sgn = GA.sign[I, J];
					r[K] += Expr.Integer(sgn) * c[I] * B.c[J];
				}
			}
		}
		return new MV(r);
	}
	// --- Grade utilities ---
	public Expr[] Grade(int grade) {
		var indices = GA.gradeMap[grade];
		var comps = new Expr[indices.Length];
		for (int i = 0; i < indices.Length; i++)
			comps[i] = c[indices[i]];
		return comps;
	}
	public MV Sym(MV B) {
		return this*B/Expr.Integer(2) + B*this/Expr.Integer(2);
	}
	public MV ASym(MV B) {
		return this*B/Expr.Integer(2) - B*this/Expr.Integer(2);
	}
	//.RationalSimplify("x")
	
	public MV Conj() {
		var r = new Expr[GA.size];
		for (int I = 0; I < GA.size; I++) {
			int g = BitOperations.PopCount((uint)I); // grade
			// reverse GA.sign = (-1)^(g*(g-1)/2)
			int sgn = ((g * (g - 1) / 2) % 2 == 0) ? +1 : -1;
			r[I] = Expr.Integer(sgn) * c[I];
		}
		return new MV(r);
	}
	public bool IsInvertible() {
		var rev = this.Conj();
		var prod = this * rev;
		// Check that all non-scalar components are zero
		for (int i = 1; i < GA.size; i++) {
			if (!prod.c[i].IsZero()) return false;  // not purely scalar → no simple inverse
		}
		// Scalar part must be non-zero
		return !prod.c[0].IsZero();
	}

	public MV Inv() {
		var rev = this.Conj();
		var prod = this * rev;
		// Check scalar-only
		for (int i = 1; i < GA.size; i++) {
			if (!prod.c[i].IsZero())
				throw new InvalidOperationException(
					"Multivector is not a versor; inverse not defined by M~ / (M M~)."
				);
		}
		var norm = prod.c[0];
		if (norm.IsZero())
			throw new InvalidOperationException("Multivector not invertible (zero norm).");

		return rev / norm;
	}

	// --- Operators ---
	public static MV operator +(MV a, MV b) {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) r[i] = (a.c[i] + b.c[i]);
		return new MV(r);
	}
	public static MV operator *(Expr s, MV a) {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) r[i] = (s * a.c[i]);
		return new MV(r);
	}
	public static MV operator -(MV a) {
		return (-1)*a;
	}
	public static MV operator -(MV a, MV b) {
		return a + (-1)*b;
	}
	public static MV operator *(MV a, Expr s) {
		return s*a;
	}
	public static MV operator /(MV a, Expr s) {
		return (Expr.Integer(1)/s)*a;
	}
	public static MV operator /(Expr s, MV a) {
		return s*a.Inv();
	}
	public static MV operator /(MV A, MV B) {
		return A*B.Inv();
	}
	
	
	public MV Simplify() {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++)
			r[i] = c[i].SimplifyRat();
		return new MV(r);
	}
	public MV SimplifyGar() {
		var r = new Expr[GA.size];
		for (int i = 0; i < GA.size; i++) {
			r[i] = c[i];
			var stra = r[i].CollectIdentifiers(); //do not pass nulls into COllectIdentifier
			foreach (var str in stra)
				r[i] = r[i].RationalSimplify(str);
		}
		return new MV(r);
	}

	
	public bool IsConstant(int n, int m, out MV t,
		double min = 0.1, double max = 2, double tol = 1e-9)
	{
		t = new MV(new Expr[GA.size]);
		for (int i = 0; i < GA.size; i++) {
			if (!c[i].IsConstant(n, m, out Expr r, min, max, tol)) {
				return false;
			}
			t.c[i] = r;
		}
		return true;
	}
	public static double[] da(int[] grades, params MV[] mvs) {
		if (grades.Length != mvs.Length)
			throw new ArgumentException("Range must have min < max.");
		int n = mvs.Length;
		double[][] res = new double[n][];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < GA.gradeMap[grades[i]].Length; j++) {
				if (mvs[i].c[GA.gradeMap[grades[i]][j]].ToDouble(out double d))
					res[i][j] = d;
			}
		}
		return res.SelectMany(a => a).ToArray();
	}
	
	

	public override string ToString() {
		return this.ToString(1);
	}
	//not working
	public string ToString(int startIndex) {
		var sb = new System.Text.StringBuilder();
		bool first = true;
		// Flatten GA.gradeMap into the dual-aware order
		var flat = GA.gradeMap.SelectMany(g => g).ToArray();
		foreach (var idx in flat) {
			var simple = c[idx];
			if (simple.IsZero()) continue;
			if (!first) sb.Append(" + \n");
			sb.Append("[");
			sb.Append(simple);
			sb.Append("]");
			if (idx != 0) {
				// always print explicit indices from GA.mask
				string mask = Convert.ToString(idx, 2).PadLeft(GA.n, '0');
				char[] chars = mask.ToCharArray();
				Array.Reverse(chars);
				var blade = new System.Text.StringBuilder();
				for (int bit = 0; bit < chars.Length; bit++) {
					if (chars[bit] == '1')
						blade.Append(bit + startIndex);
				}
				if (blade.Length > 0) {
					if (GA.bladeSign[idx] < 0) sb.Append("*-e");
					else sb.Append("*e");
					sb.Append(blade.ToString());
				}
			}

			first = false;
		}

		return first ? "0" : sb.ToString();
	}
}
