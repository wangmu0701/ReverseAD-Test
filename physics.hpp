/****
 * AD fuer QM2-RG Fluss mit Taylorloesung
 * proof of concept
 * Mathias Wagner
 * compile and link against ADOL-C
 * g++ -I $ADOLC/include/ -L$ADOLC/lib/ -ladolc cflow.cpp
 *
 */

#include <cmath>
#include <iostream>

#define PRINT_RESULTS
#define DERIVATIVE_ORDER 5
#define NUM_IND 6

/**********************************************************/
/* C language definitions for use with Mathematica output */

#define Sqrt(x)		(sqrt((x)))
#define sqr(x)		(pow((x),2))

#define Power(x, y)	(pow((x), (y)))
#define Abs(x)		(fabs((x)))

#define Exp(x)		(exp((x)))
#define Log(x)		(log((x)))

#define Sin(x)		(sin((x)))
#define Cos(x)		(cos((x)))
#define Tan(x)		(tan((x)))

#define ArcSin(x)       (asin((x)))
#define ArcCos(x)       (acos((x)))
#define ArcTan(x)       (atan((x)))

#define Sinh(x)          (sinh((x)))
#define Cosh(x)          (cosh((x)))
#define Tanh(x)          (tanh((x)))

#define Coth(x)          (1./tanh((x)))
#define Csch(x)          (1./sinh((x)))
#define Sech(x)          (1./cosh((x)))



void get_initials(double* x, int n) {
  x[0] = 9;
  x[1] = 0;
  x[2] = 50.5;
  x[3] = 0;
  x[4] = 1200; // T0;
  x[5] = 0; // mu
}

template <typename Type>
void ady(Type* yin, Type* dy, Type k, double T, Type mu)
{
	// differential equation

	// parameters, constant for now
	//double T = 10.;
	// mu is always zero but we want to calculate derivatives w.r.t. mu
	//double mu = 0.;
	double g = 3.2;
	double h = 1.77e6;
	double nu = 12;

	Type r0k = yin[0];
	//double a0k = yin[1];
	Type a4k = yin[2];
	Type a6k = yin[3];
	//double  nu;

	//Function definitions

	double Pi = 3.1415926535897932385;
	/*Coth = lambda x : 1 / math.tanh(x)
	 Tanh = math.tanh
	 Csch = lambda x: 1 / math.sinh(x)*/

	//Abkuerzungen
	Type k2 = k * k;
	Type k4 = k2 * k2;
	Type a2 = h / (sqrt(2. * r0k));
	const double g2 = g * g;
	const double g4 = g2 * g2;
	const double g6 = g4 * g2;
	const double b2 = 0.5 / T;
	Type vf = k4 * 0.008443431970194814286989955;

	const Type a46 = 3 * a4k + 2 * a6k * r0k;
	const Type a4a2r0 = (a4k + a2 / (2 * r0k));

	Type Eq2 = k * k + 2 * g2 * r0k;
	Type Eq = sqrt(Eq2);
	Type Eqm = (Eq - mu) * b2;
	Type Eqp = (Eq + mu) * b2;
        Type g2tanqq = g2 * nu * (Tanh(Eqm) + Tanh(Eqp));

	Type E12 = a2 + k2;
	Type E1 = sqrt(E12);

	const Type E22 = a2 + k2 + 2 * a4k * r0k;
	const Type E2 = sqrt(E22);

	dy[0]
			= -(vf * ((-3 * a4k * Coth(b2 * E1)) / (2. * (E12 * E1)) - (a46
					* Coth(b2 * E2)) / (2. * (E22 * E2)) - (3 * a4k
					* sqr(Csch(b2 * E1))) / (4. * E12 * T) - (a46
					* sqr(Csch(b2 * E2))) / (4. * E22 * T) - (nu * b2 * g2
					/ Eq2 * ((sqr(Sech(Eqm))) + (sqr(Sech(Eqp))))) + g2tanqq
					/ (Eq2 * Eq))) / a4a2r0;

	dy[1] = (k4 * ((3 * Coth(b2 * E1)) / E1 + Coth(b2 * E2) / E2 - (g2tanqq
			/ (g2 * Eq)))) / (12. * sqr(Pi)) - (a2 * vf * ((-3 * a4k
			* Coth(b2 * E1)) / (2. * (E12 * E1)) - (a46 * Coth(b2 * E2)) / (2.
			* (E22 * E2)) - (3 * a4k * sqr(Csch(b2 * E1))) / (4. * E12 * T)
			- (a46 * sqr(Csch(b2 * E2))) / (4. * E22 * T) - (b2 * g2 / Eq2 * nu
			* (sqr(Sech(Eqm)) + sqr(Sech(Eqp)))) + (g2tanqq) / (Eq2 * Eq)))
			/ (a4a2r0);

	dy[2] = vf * (-(a6k * ((-3 * a4k * Coth(b2 * E1)) / (2. * (E12 * E1))
			- (a46 * Coth(b2 * E2)) / (2. * (E22 * E2)) - (3 * a4k
			* sqr(Csch(b2 * E1))) / (4. * E12 * T) - (a46
			* sqr(Csch(b2 * E2))) / (4. * E22 * T) - (nu * ((b2 * g2
			* sqr(Sech(Eqm))) / Eq + (b2 * g2 * sqr(Sech(Eqp))) / Eq))
			/ Eq + g2tanqq / (Eq2 * Eq))) / (a4k + a2 / r0k) + ((3 * ((3
			* sqr(a4k)) / (4. * (E12 * E12 * E1)) - a6k / (2. * (E12 * E1)))
			* Coth(b2 * E1) + ((3 * sqr(a46)) / (4. * Power(E22, 2.5)) - (5
			* a6k) / (2. * (E22 * E2))) * Coth(b2 * E2) + (3 * sqr(a4k)
			* sqr(Csch(b2 * E1))) / (4. * sqr(E12) * T) + (3 * ((sqr(a4k)
			* sqr(Csch(b2 * E1))) / (8. * (E12 * E1) * T) - (a6k
			* sqr(Csch(b2 * E1))) / (4. * E1 * T) + (sqr(a4k)
			* Coth(b2 * E1) * sqr(Csch(b2 * E1))) / (8. * E12 * (T * T))))
			/ E1\
 + (sqr(a46) * sqr(Csch(b2 * E2))) / (4. * sqr(E22) * T)
			+ ((sqr(a46) * sqr(Csch(b2 * E2))) / (8. * (E22 * E2) * T) - (5
					* a6k * sqr(Csch(b2 * E2))) / (4. * E2 * T) + (sqr(a46)
					* Coth(b2 * E2) * sqr(Csch(b2 * E2))) / (8. * E22
					* (T * T))) / E2 - nu * ((-2 * g2 * ((b2 * g2
			* sqr(Sech(Eqm))) / Eq + (b2 * g2 * sqr(Sech(Eqp))) / Eq))
			/ (Eq2 * Eq) + (3 * g2 * g2 * (Tanh(Eqm) + Tanh(Eqp)))
			/ Power(Eq2, 2.5) - (g2 * g2 * (sqr(Sech(Eqm)) * (T + Eq
			* Tanh(Eqm)) + sqr(Sech(Eqp)) * (T + Eq * Tanh(Eqp)))) / (2.
			* (Eq2 * Eq2) * (T * T))))));

	dy[3] = vf * (3 * ((-15 * Power(a4k, 3)) / (8. * Power(E12, 3.5)) + (9
			* a4k * a6k) / (4. * Power(E12, 2.5))) * Coth(b2 * E1) + ((-15
			* Power(a46, 3)) / (8. * Power(E22, 3.5)) + (45 * a46 * a6k) / (4.
			* Power(E22, 2.5))) * Coth(b2 * E2) - (9 * a4k * ((3
			* Power(a4k, 2)) / (4. * Power(E12, 2.5)) - a6k / (2.
			* Power(E12, 1.5))) * sqr(Csch(b2 * E1))) / (4. * E1 * T) - (9
			* a4k * ((sqr(a4k) * sqr(Csch(b2 * E1))) / (8. * (E12 * E1) * T)
			- (a6k * sqr(Csch(b2 * E1))) / (4. * E1 * T) + (sqr(a4k)
			* Coth(b2 * E1) * sqr(Csch(b2 * E1))) / (8. * E12 * sqr(T))))
			/ (2. * (E12 * E1)) + (3 * ((-3 * Power(a4k, 3)
			* sqr(Csch(b2 * E1))) / (16. * Power(E12, 2.5) * T) + (3 * a4k
			* a6k * sqr(Csch(b2 * E1))) / (8. * Power(E12, 1.5) * T) - (3
			* Power(a4k, 3) * Coth(b2 * E1) * sqr(Csch(b2 * E1))) / (16.
			* sqr(E12) * sqr(T)) + (3 * a4k * a6k * Coth(b2 * E1)
			* sqr(Csch(b2 * E1))) / (8. * E12 * sqr(T)) - (Power(a4k, 3)
			* sqr(Coth(b2 * E1)) * sqr(Csch(b2 * E1))) / (16. * (E12 * E1)
			* Power(T, 3)) - (Power(a4k, 3) * Power(Csch(b2 * E1), 4)) / (32.
			* (E12 * E1) * Power(T, 3)))) / E1 - (3 * a46 * ((3 * sqr(a46))
			/ (4. * Power(E22, 2.5)) - (5 * a6k) / (2. * (E22 * E2)))
			* sqr(Csch(b2 * E2))) / (4. * E2 * T) - (3 * a46 * ((sqr(a46)
			* sqr(Csch(b2 * E2))) / (8. * (E22 * E2) * T) - (5 * a6k
			* sqr(Csch(b2 * E2))) / (4. * E2 * T) + (sqr(a46)
			* Coth(b2 * E2) * sqr(Csch(b2 * E2))) / (8. * E22 * sqr(T))))
			/ (2. * E22 * E2) + ((-3 * Power(a46, 3) * sqr(Csch(b2 * E2)))
			/ (16. * Power(E22, 2.5) * T) + (15 * a46 * a6k
			* sqr(Csch(b2 * E2))) / (8. * (E22 * E2) * T) - (3
			* Power(a46, 3) * Coth(b2 * E2) * sqr(Csch(b2 * E2))) / (16.
			* sqr(E22) * sqr(T)) + (15 * a46 * a6k * Coth(b2 * E2)
			* sqr(Csch(b2 * E2))) / (8. * E22 * sqr(T)) - (Power(a46, 3)
			* sqr(Coth(b2 * E2)) * sqr(Csch(b2 * E2))) / (16. * (E22 * E2)
			* Power(T, 3)) - (Power(a46, 3) * Power(Csch(b2 * E2), 4)) / (32.
			* Power(E22, 1.5) * Power(T, 3))) / E2 - nu * ((9 * g4 * ((b2 * g2
			* sqr(Sech(Eqm))) / Eq + (b2 * g2 * sqr(Sech(Eqp))) / Eq))
			/ Power(Eq2, 2.5) - (15 * g6 * (Tanh(Eqm) + Tanh(Eqp)))
			/ Power(Eq2, 3.5) - (3 * g2 * (-(g4 * sqr(Sech(Eqm))) / (2.
			* Power(Eq2, 1.5) * T) - (g4 * sqr(Sech(Eqp))) / (2.
			* Power(Eq2, 1.5) * T) - (g4 * sqr(Sech(Eqm)) * Tanh(Eqm)) / (2.
			* Eq2 * Power(T, 2)) - (g4 * sqr(Sech(Eqp)) * Tanh(Eqp)) / (2.
			* Eq2 * Power(T, 2)))) / Power(Eq2, 1.5)\
 + ((3 * g6
			* sqr(Sech(Eqm))) / (2. * Power(Eq2, 2.5) * T) - (g6
			* Power(Sech(Eqm), 4)) / (4. * Power(Eq2, 1.5) * Power(T, 3))
			+ (3 * g6 * sqr(Sech(Eqp))) / (2. * Power(Eq2, 2.5) * T) - (g6
			* Power(Sech(Eqp), 4)) / (4. * Power(Eq2, 1.5) * Power(T, 3))
			+ (3 * g6 * sqr(Sech(Eqm)) * Tanh(Eqm)) / (2. * Power(Eq2, 2)
					* sqr(T)) + (g6 * sqr(Sech(Eqm)) * sqr(Tanh(Eqm))) / (2.
			* (Eq2 * Eq) * Power(T, 3)) + (3 * g6 * sqr(Sech(Eqp)) * Tanh(Eqp))
			/ (2. * sqr(Eq2) * sqr(T)) + (g6 * sqr(Sech(Eqp)) * sqr(Tanh(Eqp)))
			/ (2. * (Eq2 * Eq) * Power(T, 3))) / Eq));

}

template <typename T>
T eval_func(T* xad, int n) {
  T yad[4];
  ady(xad, yad, xad[5], 150, xad[6]);
  T fad = yad[0] + yad[1] + yad[2] + yad[3];
  return fad;
}
