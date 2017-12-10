#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using namespace std;
using namespace Eigen;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = (j[1]["speed"]);
          double steer_value= (j[1]["steering_angle"]);
          double throttle_value= (j[1]["throttle"]);
// From miles per hour to meters per second
v=v*0.44704;

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].

          //Display the MPC predicted trajectory 

		vector<double> mpc_x_vals=ptsx;
		vector<double> mpc_y_vals=ptsy;
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

       //   msgJson["mpc_x"] = mpc_x_vals;
       //   msgJson["mpc_y"] = mpc_y_vals;
vector<double> waypoints_car_x(ptsx.size());
vector<double> waypoints_car_y(ptsx.size());
vector<double> dif_x(ptsx.size());
vector<double> dif_y(ptsx.size());
VectorXd wayx(ptsx.size());
VectorXd wayy(ptsx.size());
//cout << "hei hei"<< endl;
for(unsigned int i=0;i<ptsx.size();i++){
dif_x[i]=(ptsx[i]-px);
dif_y[i]=(ptsy[i]-py);
}

for(unsigned int i=0;i<ptsx.size();i++){
waypoints_car_x[i]=dif_x[i]*cos(-psi)-dif_y[i]*sin(-psi);
waypoints_car_y[i]=dif_y[i]*cos(-psi)+dif_x[i]*sin(-psi);

// Previously I had
//waypoints_car_x[i]=dif_x[i]*cos(psi)-dif_y[i]*sin(-psi);
//waypoints_car_y[i]=dif_y[i]*cos(psi)+dif_x[i]*sin(-psi);

wayx[i]=waypoints_car_x[i];
wayy[i]=waypoints_car_y[i];
}

  Eigen::VectorXd coefs2 =polyfit(wayx,wayy,3);
  cout << "coefs " << coefs2 << endl;

vector<double> wayx_intr;
vector<double> wayy_intr;

//for(int i=0;i<20;i=i+4){
//wayx_intr.push_back(double(i));
//wayy_intr.push_back(polyeval(coefs2,i));
//}
double cte=coefs2[0];
double epsi= -atan(coefs2[1]);

//double 
//steer_value=msgJson["steering_angle"];
//double 
//throttle_value=msgJson["throttle"]; 
const double Lf = 2.67;
int latency_ms=100;
//int latency_ms=100;
double dt=latency_ms/1000;
//double LF=Lf;
double x_new=v*dt;
double y_new=0;
double psii_new=-v/Lf*steer_value*dt;

double v_new=v+throttle_value*dt;
double cte_new = cte+v*sin(epsi)*dt;
double epsi_new =epsi + psii_new;// v/Lf*steer_value*dt;
//steer_value=-cte*0.3;
//throttle_value=0.4;

VectorXd state_after_latency(6);
state_after_latency << x_new,y_new,psii_new,v_new,cte_new,epsi_new;

vector <double> mpc_optimized=mpc.Solve(state_after_latency,coefs2);

// Change then also from the MPC equations
steer_value=mpc_optimized[0]/deg2rad(25);
//steer_value=-mpc_optimized[0]/deg2rad(25);


throttle_value=mpc_optimized[1];
cout << "Steering result: steer" << steer_value << "throttle_value" << throttle_value <<endl;

mpc_x_vals=mpc.plottable_trajectory_x;
mpc_y_vals=mpc.plottable_trajectory_y;
msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;
msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;


vector<double> next_x_vals=wayx_intr;
vector<double> next_y_vals=wayy_intr;


          //Display the waypoints/reference line     
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

 //         msgJson["next_x"] = next_x_vals;
  //        msgJson["next_y"] = next_y_vals;
          msgJson["next_x"] = waypoints_car_x;
          msgJson["next_y"] = waypoints_car_y;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
      
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
           this_thread::sleep_for(chrono::milliseconds(latency_ms));
//this_thread::sleep_for(chrono::milliseconds(0));

          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
