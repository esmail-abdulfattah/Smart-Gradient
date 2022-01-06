
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <LBFGS.h>
#include <random>

std::normal_distribution<double> distribution(0.0,7e-8);
std::default_random_engine generator;
double mach_eps = std::numeric_limits<double>::epsilon();

using namespace LBFGSpp;
using namespace Eigen;
using Eigen::MatrixXd;
typedef double (* fun)(const VectorXd& x); 

class MyFunction
{
    private:
        Eigen::MatrixXd G, c_G;
        Eigen::VectorXd prev_x, curr_x;
        int n, count;
    public:
        MyFunction(int n_) : n(n_) {
            G.resize(n,n); c_G.resize(n,n);
            for(size_t i=0;i<n;i++) {G(i,i) = 1.0; c_G(i,i) = 1.0;}
            prev_x = VectorXd::Zero(n);
            curr_x = VectorXd::Zero(n);
            count = 0;
        }
        double myfun(const VectorXd& x)
        {
            int n = x.size();
            double fx = 0.0;
            for(int i = 0; i < n; i += 2)
            {
                double t1 = 1.0 - x[i];
                double t2 = 10 * (x[i + 1] - x[i] * x[i]);
                
                fx += t1 * t1 + t2 * t2;
            }
            return fx;
        }
        double mytfun(const VectorXd& phi)
        {
            VectorXd xx = curr_x + G*phi;
            return myfun(xx);
        }
        void mytgrad(VectorXd& grad){
            int n = grad.size();

            double step = 0.0001;
            for (size_t i = 0; i < n; i++)
            {
                VectorXd x = VectorXd::Zero(n); //gradient always should be evaluated at zero
                double y1=0.0,y2=0.0;
                x[i] += step;
                y1=mytfun(x);
                x[i] -= 2*step;
                y2=mytfun(x);
                grad[i] = (y1 - y2) / (2.0*step);
            }}

            
        void scale(VectorXd &x){  
            double m = x.mean() ;
            double sd = 0.0;
            for(size_t i=0;i<x.size();i++) sd += (x[i] - m)*(x[i] - m);
            sd = sqrt(sd/(x.size()-1.0));
            for(size_t i=0;i<x.size();i++)  x[i] = (x[i] - m)/sd;}

        void set_prev_x(VectorXd x)  {prev_x = x; count++;}
        void update_G(VectorXd &current_x)  {
            curr_x = current_x;
            for(size_t i=(n-1);i>=1;i--) for(size_t j=0;j<n;j++) c_G(j,i) = c_G(j,i-1);
            VectorXd xdiff = current_x - prev_x ; 
            for(size_t ii=0; ii<n; ii++) xdiff[ii] = xdiff[ii] + distribution(generator);
            scale(xdiff);
            for(size_t i=0;i<n;i++) c_G(i,0) = xdiff[i] ;
            G = c_G;}  
        void MGS_orthogonalization(){   
            size_t i, j, k;
            double r;
            VectorXd q = VectorXd::Zero(n);
            for (i = 0; i < n; i++) 
            {
                for (r = 0.0, j = 0; j < n; j++) r += (G(j,i)*G(j,i));
                r = std::sqrt(r);
                for (j = 0; j < n; j++) {q[j] = G(j,i)/r; G(j,i) = q[j];}

                for (j = i + 1; j < n; j++) 
                {
                    for (r = 0, k = 0; k < n; k++) r += q[k] * G(k,j); 
                    for (k = 0; k < n; k++) G(k,j) = G(k,j) - r*q[k]; 
                }
            }
        }
        void get_grad(VectorXd &grad) {grad = (G.transpose()).colPivHouseholderQr().solve(grad);}
        void transform_grad(VectorXd &x){
            if(count>0) update_G(x); else G = c_G;
            MGS_orthogonalization();
            set_prev_x(x);}

        void smartGrad(const VectorXd& x,VectorXd& grad){
            VectorXd xx = x;
            transform_grad(xx);
            mytgrad(grad);
            get_grad(grad);}
        void vanillaGrad(const VectorXd& xx,VectorXd& grad){
            int n = grad.size();

            double step = 0.0001;
            for (size_t i = 0; i < n; i++)
            {
                VectorXd x = xx; 
                double y1=0.0,y2=0.0;
                x[i] += step;
                y1=myfun(x);
                x[i] -= 2*step;
                y2=myfun(x);
                grad[i] = (y1 - y2) / (2.0*step);
            }
        }
        void exactgradient(const VectorXd& x, VectorXd& grad)
        {
            for(int i = 0; i < n; i += 2)
            {
                double t1 = 1.0 - x[i];
                double t2 = 10 * (x[i + 1] - x[i] * x[i]);
                grad[i + 1] = 20 * t2;
                grad[i] = -2.0 * (x[i] * grad[i + 1] + t1);
            }
        }
        double operator()(const VectorXd& x, VectorXd& grad)
        {
            smartGrad(x,grad);
            return myfun(x);
        }

};

int main()
{
    const int n = 4;
    LBFGSParam<double> param;
    //param.epsilon = 0.001;
    param.max_iterations = 100;
    LBFGSSolver<double,LineSearchBacktracking> solver(param);
    MyFunction fun(n);

    VectorXd x = VectorXd::Zero(n);
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;


    return 0;
}

