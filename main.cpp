#include <iostream>
#include <vector>
#include <clocale>


using namespace std;

struct polynomial {
    vector<int> monomials;
    int degree = int(monomials.size()) - 1;
};

struct success_and_result {
    bool success = false;
    polynomial result;
};

#include "pfadd.cpp"

int calculate(const polynomial &f, int n) {
    int sum = 0;
    int p = 0;
    for (auto i : f.monomials) {
        sum += i * int(pow(n, p++));
    }
    return sum;
}

bool equal(double x, double y) {
    return fabs(x - y) < 1e-7;
}

success_and_result divide(polynomial A, polynomial B) {
    int n = (int) A.monomials.size() - 1;
    int m = (int) B.monomials.size() - 1;
    polynomial Q{vector<int>(n - m + 1)};
    for (int i = n; i >= m; i--) {
        Q.monomials[i - m] = A.monomials[i] / B.monomials[m];
        for (int j = m; j >= 0; j--)
            A.monomials[i - m + j] -= B.monomials[j] * Q.monomials[i - m];
    }
    while (!A.monomials.empty() && equal(A.monomials.back(), 0))
        A.monomials.pop_back();
    if (A.monomials.empty())
        return success_and_result{true, Q};
    else
        return success_and_result{false};
}

vector<int> find_divisors(int n) {
    vector<int> ans;
    int i;
    ans.push_back(1);
    ans.push_back(-1);
    for (i = 2; i * i < n; i++) {
        if (n % i == 0) {
            ans.push_back(i);
            ans.push_back(-i);
            ans.push_back(n / i);
            ans.push_back(-n / i);
        }
    }
    if (i * i == n) {
        ans.push_back(i);
        ans.push_back(-i);
    }
    ans.push_back(n);
    ans.push_back(-n);
    return ans;
}

vector<vector<int>> direct_multiplication(const vector<vector<int>> &A, const vector<int> &B) {
    vector<vector<int>> ans;
    for (const auto &a : A) {
        for (const int &b : B) {
            auto t(a);
            t.push_back(b);
            ans.push_back(t);
        }
    }
    return ans;
}

vector<double> get_lagrange(const vector<double> &x, vector<double> f) {
    int N = x.size();
    vector<double> c(N), temp(N);

    c[0] = f[0];
    for (int i = 1; i < N; i++) {
        for (int j = 0; j < N - i; j++) temp[j] = (f[j + 1] - f[j]) / (x[j + i] - x[j]);
        f = temp;
        c[i] = f[0];
    }
    return c;
}

vector<double> standard_polynomial(const vector<double> &c, const vector<double> &x) {
    int N = x.size();
    vector<double> a(N, 0.0);
    vector<double> p(N), prev(N);

    p[0] = 1;
    a[0] = c[0] * p[0];
    for (int i = 1; i < N; i++) {
        prev = p;
        p[0] = -x[i - 1] * prev[0];
        a[0] += c[i] * p[0];
        for (int j = 1; j <= i; j++) {
            p[j] = prev[j - 1] - x[i - 1] * prev[j];
            a[j] += c[i] * p[j];
        }
    }

    return a;
}

success_and_result interpol_lagrange(const vector<double> &x, const vector<double> &y) {
    vector<double> c = get_lagrange(x, y);
    vector<double> arrayp = standard_polynomial(c, x);

    polynomial ans;
    for (double i : arrayp) {
        if (!equal(i, int(i))) {
            return success_and_result{false};
        } else {
            ans.monomials.push_back(int(i));
        }
    }
    return success_and_result{true, ans};
}

void print_polynomial(polynomial a) {
    int n = a.degree;
    if (n > 1) {
        if (a.monomials[n] == -1)
            cout << "-x^" << n;
        else if (a.monomials[n] == 1)
            cout << "x^" << n;
        else
            cout << a.monomials[n] << "x^" << n;
    }

    for (int i = n - 1; i > 1; --i) {
        if (a.monomials[i] != 0) {
            if (a.monomials[i] > 0) {
                cout << '+';
                (a.monomials[i] == 1) ? cout << "x^" << i : cout << a.monomials[i] << "x^" << i;
            }
            else
                (a.monomials[i] == -1) ? cout << "-x^" << i : cout << a.monomials[i] << "x^" << i;
        }
    }

    if (a.monomials[1] != 0) {
        if (a.monomials[1] > 0) {
            if (n > 1)
                cout << '+';
            (a.monomials[1] == 1) ? cout << "x" : cout << a.monomials[1] << "x";
        }
        else
            (a.monomials[1] == -1) ? cout << "-x" : cout << a.monomials[1] << "x";
    }

    if (a.monomials[0] != 0) {
        if (a.monomials[0] > 0) {
            cout << '+';
            cout << a.monomials[0];
        }
        else
            cout << a.monomials[0];
    }
    cout << endl;
}

success_and_result factorize(const polynomial &f) {
    polynomial g;
    vector<int> M;
    vector<vector<int>> U;
    success_and_result interpol;


    for (int i = 0; i * 2 <= f.degree; i++) {
        if (calculate(f, i) == 0) {
            g.monomials = vector<int>{-i, 1};
            g.degree = 1;
            return success_and_result{true, g};
        }
    }
    auto f0 = calculate(f, 0);
    auto divisors_f0 = find_divisors(f0);
    U = direct_multiplication(vector<vector<int>>(1), divisors_f0);
    for (int i = 1; i * 2 <= f.degree; i++) {
        auto fi = calculate(f, i);
        M = find_divisors(fi);
        U = direct_multiplication(U, M);
        for (const vector<int> &u : U) {
            auto poly_u = polynomial{u};

            vector<double> x, y;

            for (int j = 0; j <= i; j++) {
                x.push_back(double(j));
                y.push_back(double(calculate(poly_u, j)));
            }

            interpol = interpol_lagrange(x, y);
            // print_vec<double>(x);
            // print_vec<double>(y);
            // cout << "success = " << interpol.success << '\n';
            // print_polynomial(poly_u);

            if (interpol.success) {
                g = interpol.result;
                // print_polynomial(g);
                if (divide(f, g).success) {
                    g.degree = i;
                    return success_and_result{true, g};
                }
            }
        }
    }
    return success_and_result{false};
}


/*
int main()
{
    vector<double> x = { 1, 2, 3, 4 };
    vector<double> y = { 6, 9, 2, 5 };
    //vector<vector<double>> g_points(2);
    //g_points[0] = x;
    //g_points[1] = y;
    auto ip = interpol_lagrange(x, y);

//    cout << "Newton polynomial:   ";   writeNewtonPolynomial( c, x );
    cout << "Standard polynomial: ";   print_polynomial( ip.result );

}
*/
int main() {
    int n;
    setlocale(LC_ALL, "Russian");
    cout << "¬ведите степень многочлена:\n";
    cin >> n;
    cout << "¬ведите одночлены (от наибольшей степени к наименьшей):\n";
    polynomial a{vector<int>(n + 1), n};
    for (int i = n; i >= 0; --i) {
        cin >> a.monomials[i];
    }

    cout << "¬веденный многочлен:\n";
    print_polynomial(a);

    auto factorized = factorize(a);
    vector<polynomial> answer;
    while (factorized.success) {
        answer.push_back(factorized.result);
        a = divide(a, factorized.result).result;
        factorized = factorize(a);
    }

    if (!a.monomials.empty()) {
        answer.push_back(a);
    }

    cout << "ќтвет:\n";
    for (const auto& poly : answer) {
        print_polynomial(poly);
    }
    return 0;
}


