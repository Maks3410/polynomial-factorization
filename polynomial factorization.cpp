//#include <iostream>
//#include <vector>
//
//using namespace std;
//
//struct polynomial {
//    vector<int> monomials;
//    int degree = int(monomials.size()) - 1;
//};
//
//struct success_and_result {
//    bool success = false;
//    polynomial result;
//};
//
//int calculate(const polynomial& f, int n) {
//    int sum = 0;
//    int p = 0;
//    for (auto i : f.monomials) {
//        sum += i * int(pow(n, p++));
//    }
//    return sum;
//}
//
//bool equal(double x, double y) {
//    return fabs(x - y) < 1e-9;
//}
//
//success_and_result divide(polynomial A, polynomial B) {
//    int n = (int)A.monomials.size() - 1;
//    int m = (int)B.monomials.size() - 1;
//    polynomial Q{vector<int>(n - m + 1)};
//    for (int i = n; i >= m; i--) {
//        Q.monomials[i - m] = A.monomials[i] / B.monomials[m];
//        for (int j = m; j >= 0; j--)
//            A.monomials[i - m + j] -= B.monomials[j] * Q.monomials[i - m];
//    }
//    A.monomials.resize(m);
//    while(!A.monomials.empty() && equal(A.monomials.back(), 0))
//        A.monomials.pop_back();
//    if (A.monomials.empty())
//        return success_and_result{true, Q};
//    else
//        return success_and_result{false};
//}
//
//vector<int> find_divisors(int n) {
//    vector<int> ans;
//    int i;
//    ans.push_back(1);
//    ans.push_back(-1);
//    for (i = 2; i * i < n; i++) {
//        if (n % i == 0) {
//            ans.push_back(i);
//            ans.push_back(-i);
//            ans.push_back(n / i);
//            ans.push_back(-n / i);
//        }
//    }
//    if (i * i == n) {
//        ans.push_back(i);
//        ans.push_back(-i);
//    }
//    ans.push_back(n);
//    ans.push_back(-n);
//    return ans;
//}
//
//vector<vector<int>> direct_multiplication(const vector<vector<int>>& A, const vector<int>& B) {
//    vector<vector<int>> ans;
//    for (const auto &a : A) {
//        for (const int & b : B) {
//            auto t(a);
//            t.push_back(b);
//            ans.push_back(t);
//        }
//    }
//    return ans;
//}
//
//vector<double> getNewton( const vector<double> &x, vector<double> f )          // Gets Newton polynomial
//{
//    int N = x.size();
//    vector<double> c(N), temp(N);
//
//    c[0] = f[0];
//    for ( int i = 1; i < N; i++ )       // Compute ith differences
//    {
//        for ( int j = 0; j < N - i; j++ ) temp[j] = ( f[j+1] - f[j] ) / ( x[j+i] - x[j] );
//        f = temp;
//        c[i] = f[0];
//    }
//    return c;
//}
//
//vector<double> standardPolynomial( const vector<double> &c, const vector<double> &x )    // Multiplies out Newton polynomial
//{
//    int N = x.size();
//    vector<double> a( N, 0.0 );                   // a[j] holds coefficient of x^j in final sum
//    vector<double> p(N), prev(N);                 // At the ith step, p[ ] is the ith polynomial of Newton factors
//
//    p[0] = 1;
//    a[0] = c[0] * p[0];
//    for ( int i = 1; i < N; i++ )
//    {
//        prev = p;
//        p[0] = -x[i-1] * prev[0];
//        a[0] += c[i] * p[0];
//        for ( int j = 1; j <= i; j++ )
//        {
//            p[j] = prev[j-1] - x[i-1] * prev[j];
//            a[j] += c[i] * p[j];
//        }
//    }
//
//    return a;
//}
//
//success_and_result interpol_newton(vector<vector<double>> Net){
//    vector<double> c = getNewton( Net[0], Net[1] );
//    vector<double> arrayp = standardPolynomial( c, Net[0] );
//
//    polynomial ans;
//    bool f = true;
//    for (int i = 0; i < arrayp.size(); ++i) {
//        if (!equal(arrayp[i], int(arrayp[i]))) {
//            f = false;
//        }
//        else {
//            ans.monomials.push_back(int(arrayp[i]));
//        }
//    }
//    return success_and_result{f, ans};
//}
//
//success_and_result factorize(const polynomial& f) {
//    polynomial g;
//    vector<int> M;
//    vector<vector<int>> U;
//    // vector<vector<double>> g_points(2);
//    success_and_result interpol;
//
//    bool success = false;
//
//    for (int i = 0; i * 2 <= f.degree && !success; i++) {
//        if (calculate(f, i) == 0) {
//            success = true;
//            g.monomials = vector<int>{-i, 1};
//            g.degree = 1;
//        }
//    }
//    if (!success) {
//        auto f0 = calculate(f, 0);
//        auto divisors_f0 = find_divisors(f0);
//        U = direct_multiplication(vector<vector<int>>(1), divisors_f0);
//        for(int i = 1; i * 2 <= f.degree; i++) {
//            auto fi = calculate(f, i);
//            M = find_divisors(fi);
//            U = direct_multiplication(U, M);
//            for (const vector<int>& u : U) {
//                auto poly_u = polynomial{u};
//                vector<vector<double>> g_points(2);
//                for (int j = 0; j <= i; j++) {
//                    g_points[0].push_back(double(j));
//                    g_points[1].push_back(double(calculate(poly_u, j)));
//                }
//
//                interpol = interpol_newton(g_points);
//
//                if (interpol.success) {
//                    g = interpol.result;
//                    if (divide(f, g).success) {
//                        g.degree = i;
//                        return success_and_result{true, g};
//                    }
//                }
//            }
//        }
//    }
//    return success_and_result{false};
//}
//
//
//int main() {
//    int n;
//    cout << "Enter polynomial degree:";
//    cin >> n;
//    cout << "Enter monomials (high to low):";
//    polynomial a{vector<int>(n + 1), n};
//    for (int i = n; i >= 0; --i) {
//        cin >> a.monomials[i];
//    }
//    cout << "Answer:\n";
//
//    auto factorized = factorize(a);
//    vector<polynomial> answer;
//    while (factorized.success) {
//        answer.push_back(factorized.result);
//        a = divide(a, factorized.result).result;
//        factorized = factorize(a);
//    }
//
//    if (!a.monomials.empty()) {
//        answer.push_back(a);
//    }
//
//    for (auto polynom : answer) {
//        for (auto mono : polynom.monomials) {
//            cout << mono << ' ';
//        }
//        cout << endl;
//    }
//    cin >> n;
//    return 0;
//}
//