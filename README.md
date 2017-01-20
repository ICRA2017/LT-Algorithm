LT ALGORITHM MATLAB IMPLEMENTATION

Overview:
This set of matlab scripts are an implementation of the LT Algorithm for air flow 
modeling introduced in the article "Probabilistic Air Flow Modelling Using Turbulent 
and Laminar Characteristics for Ground and Aerial Robots", presented by Hernandez 
and co-authors. For mobile robots that operate in complex, uncontrolled environments, 
estimating air flow models can be of great importance. Aerial robots use air flow 
models to plan optimal navigation paths and to avoid turbulence-ridden areas. 
Search and rescue platforms use air flow models to infer the location of gas leaks. 
Environmental monitoring robots enrich pollution distribution maps by integrating 
the information conveyed by an air flow model.

The LT Algorithm uses wind data collected at a sparse number of locations to estimate 
joint probability distributions over wind speed and direction at given query locations. 
The LT algorithm uses a novel extrapolation approach that models the air flow as a 
linear combination of laminar and turbulent components. Experimental validation of 
the LT algorithm has shown that has a high degree of stability with respect to 
parameter selection while outperforming conventional extrapolation approaches. 

Please refer to matlab scripts "example_01_pdf_estimation" and 
"example_02_airflow_modeling" for examples on how to use the code in this repository.

When using this code, please cite the work:

V. Hernandez Bennetts, T. P. Kucner, E. Schaffernicht, P. P. Neumann, H. Fan and 
A. J. Lilienthal. Probabilistic Air Flow Modelling Using Turbulent and Laminar 
Characteristics for Ground and Aerial Robots". Robotics and Automation Letters 
(RAL). 2017. Accepted for publication.

Author contact information:
Victor Hernandez Bennetts
victor.hernandez@oru.se
http://mrolab.eu/vrbs







