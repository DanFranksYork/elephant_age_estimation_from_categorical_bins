data {
    int N; // number of individuals
    int K; // number of age categories
    int<lower=1, upper=K> age_category_index[N]; // age categories
}

parameters {
  vector<lower=0>[N] observed_age_std;
  vector<lower=0>[N] true_age;
  // real<lower=0> sigma_age[N];
  real<lower=0> sigma_age;
  real<lower=0> shape_std;
  real<lower=0> scale_std;
}

transformed parameters {
  ordered[K-1] thresholds;
  vector<lower=0>[N] observed_age;
  real<lower=0> shape;
  real<lower=0> scale;
  // Thresholds for age classes
  thresholds[1] = 15;
  thresholds[2] = 30;
  thresholds[3] = 45;
  thresholds[4] = 60;
  thresholds[5] = 75;
  // Non-centred age. The same as observed_age ~ normal(true_age,sigma_age)
  observed_age = true_age + sigma_age*observed_age_std; // if one sigma for all
  
  // In a loop to have a different sigma for each individual
  //for(i in 1:N) {
  //  observed_age[i] = true_age[i] + sigma_age[i]*observed_age_std[i];
  //}
  // Non-centred shape. The same as shape = normal(1.2,0.1)
  shape = 1.2 + 0.1*shape_std;
  // Non-centred scale. The same as scale = normal(30,1);
  scale = 30 + scale_std;
}

model {
  for(i in 1:N) {
    age_category_index[i] ~ ordered_logistic(observed_age[i], thresholds);
  }
  observed_age_std ~ std_normal();
  sigma_age ~ exponential(0.5);
  true_age ~ weibull(shape, scale);
  shape_std ~ std_normal();
  scale_std ~ std_normal();
}




//functions {
  // 1 - gompertz to give survival to an age
  // a is the asymtote
  // b shifts x axis
  // c + number, 'growth' rate
  // t number of years (first parameter is response)
//  real gompertz_lpdf(vector t, real a, real b, real c) {
//    return sum(log(1-(a * exp(-b*exp(-c*t)))));
//  }
//}
