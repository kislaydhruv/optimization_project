#include<bits/stdc++.h>
using namespace std ;


// to get the question number and marking it as a constant 
const int QUE = [](){
    int val ;
    //taking question number as input 
    cout<<"Please assign the question number "<<endl;
    cin>>val ;
    return val ;
}();




//to calculate the vlue of functions and function evaluation
double fun_eval(double& x , int& total_func_eval){
    // to increase the number of function evaluation every time it get's called 
    total_func_eval++;
    double f = 0.0 ;
    // using switch cases to assign function according to question number
    switch(QUE){
        case 1:
        f = (pow(2*x-5, 4) - pow(x*x - 1, 3))*(-1) ;
        break;
        case 2:
        f = (8+x*x*x-2*x-2*exp(x))*(-1) ;
        break;
        case 3:
        f = (-4)*x*(sin(x));
        break;
        case 4:
        f = 2*(pow(x-3, 2)) + exp(0.5*x*x);
        break;
        case 5:
        f = x*x - 10*exp(0.1*x) ;
        break;
        case 6:
        f = (20*sin(x) - 15*x*x)*(-1);
        break;
    }
    return f ;
}



//to calculate the derivative of the function numerically (central difference)
double fun_der(double&x , double& d , int& total_func_eval){
    double x_plus_d = x+d;
    double fd_a = fun_eval(x_plus_d,total_func_eval);
    double x_minus_d = x-d;
    double fd_b = fun_eval(x_minus_d,total_func_eval);
    return ((fd_a-fd_b)/(2*d)) ;
}





// Bounding Phase Method (bracketing method)
vector<double> bounding_phase(double& a , double& b , double d , double& initial_guess , int& total_func_eval){
    // we are basically doing step 1 
    double x0 = initial_guess ;
    double x_plus  = x0 + fabs(d) ;
    double x_minus = x0 - fabs(d) ;
    double f_x0    = fun_eval(x0,total_func_eval);
    double f_plus  = fun_eval(x_plus,total_func_eval);
    double f_minus = fun_eval(x_minus,total_func_eval);

    // Step 2 of the algorithm: decide the sign of the step.
    if(f_minus>=f_x0 && f_x0>=f_plus){
        d = fabs(d) ;              // f decreases towards +x -> search positive
    }
    else if(f_minus<=f_x0 && f_x0<=f_plus){
        d = -1*fabs(d) ;           // f decreases towards -x -> search negative
    }
    else{
        double j = total_func_eval ;

        return {max(a,x_minus), min(b,x_plus), j} ;
    }

    int k = 0 ;
    double x_km1 = x0 ;                 // x(k-1)
    double x_k   = x0 , f_k = f_x0 ;     // x(k)
    double x_kp1 = x_k + pow(2,k)*d ;    // x(k+1)
    double f_kp1 = fun_eval(x_kp1, total_func_eval) ;

    while(f_kp1 < f_k){
        x_km1 = x_k ;
        x_k   = x_kp1 ;
        f_k   = f_kp1 ;

        k++ ;
        x_kp1 = x_k + pow(2,k)*d ;
        f_kp1 = fun_eval(x_kp1, total_func_eval) ;
    }
    //using max , min to make sure we are returning lower bound first then the upper bound 

    double lower = min(x_km1, x_kp1) ;
    double upper = max(x_km1, x_kp1) ;
    //we are typecasting here because the return type is vector<double>
    double j = total_func_eval ;
    return {lower, upper, j} ;
}






// Bisection Method 
vector<double> bisection_method(double& a , double& b , double& d , double& e , int& total_func_eval ){
    //fingding the derivative
    double der_a = fun_der(a,d,total_func_eval);
    double der_b = fun_der(b,d,total_func_eval);


    // to check if we get some saddle point 
    if(der_a == 0.0 && der_b > 0){
        double f_a = fun_eval(a, total_func_eval);
        return {a, f_a, (double)total_func_eval};
    }
    if(der_b == 0.0 && der_a < 0){
        double f_b = fun_eval(b, total_func_eval);
        return {b, f_b, (double)total_func_eval};
    }

    // throwing error is the required criteria is not fulfilling 
    if(!(der_a < 0 && der_b > 0)){
        throw std::invalid_argument("bracket does not satisfy f'(a)<0, f'(b)>0");
    }
   
    // step 1 of bisection method 
    double z = (a+b)/2.0 ;
    double der_z = fun_der(z,d,total_func_eval);

   // subsequent steps till the algorithm get's terminated 
    while(fabs(der_z) > e){
        if(der_z > 0){
            b = z ;
        }
        else{
            a = z ;
        }
        z = (a+b)/2.0 ;
        der_z = fun_der(z,d,total_func_eval);
    }

    double f_z = fun_eval(z, total_func_eval) ;
    double j = total_func_eval;
    return {z, f_z, j};
}




int main(){
    cout<<"The question number is "<<QUE<<endl ;


    ifstream fin("input.txt");
    if(!fin){
        cerr << "Could not open input.txt (expected in the working directory)." << endl ;
        return 1 ;
    }

    double a=0,b=0,delta=0,eps=0,x0=0,h=0 ;
    bool found = false ;
    string line ;
    while(getline(fin, line)){
        if(line.empty() || line[0]=='#') continue ;
        stringstream ss(line);
        int q ; double A,B,D,E,X0,H ;
        if(ss >> q >> A >> B >> D >> E >> X0 >> H){
            if(q == QUE){
                a=A; b=B; delta=D; eps=E; x0=X0; h=H;
                found = true ;
                break ;
            }
        }
    }
    fin.close();

    if(!found){
        cerr << "No row for question " << QUE << " found in input.txt" << endl ;
        return 1 ;
    }

    int nEval = 0 ;

    //  bracket the minimum 
    vector<double> bracket = bounding_phase(a, b, delta, x0, nEval);
    double a1 = bracket[0], b1 = bracket[1];

    cout << fixed << setprecision(6);
    cout << "Bracket from bounding phase : [" << a1 << ", " << b1
         << "]  (nEval so far = " << (int)bracket[2] << ")" << endl ;

    // refine with the bisection method 
    vector<double> result = bisection_method(a1, b1, h, eps, nEval);

    cout << "x(0)   = " << x0 << endl ;
    cout << "x*     = " << result[0] << endl ;
    cout << "f(x*)  = " << result[1] << endl ;
    cout << "nEval  = " << (int)result[2] << endl ;

    return 0 ;
}
