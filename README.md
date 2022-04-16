# Rocket Propulsion Project
# Summary
This code was the final project in my rocket propulsion class at UAH. I worked with onother team member on other aspects of the project but I wrote all of this code. The project is for simulation the perforamce of solid rocket, mainly the maximum height reached of the rocket design. It calculates interal ballistics of the motor as well as atmospheric conditions through the flight which affect both the motor performace and the drag on the vehicle. After the burn phase the simulation will continue the anlyasis into the coast phase to find the apogee. The goal of this project was to create 3 different vehicle designs for 3 different altitueds. Each of these altidues needs to be achived only by changing the number of propellent grains, nozzel dimensions, and a small amount of balast. Other parameters of the vehicle we required to remain the same, this mimicked the real world need to be able to use common parts for diffrent missions to reduce cost. 

# Vehicle Design Modifications
If the parameters of the need to be changes for a diffrent mission all the parameters can be changed in src/PR08D_V2.m in the section shown below.

![image](https://user-images.githubusercontent.com/49332395/163659217-fa18c762-0680-4b3b-82cb-ed0ed8c703b9.png)

# Outputs
Many different outputs will be given for each design. There are both tables and graphs as outputs. An example of the output tables and graphs are shown below

![image](https://user-images.githubusercontent.com/49332395/163659319-51cca7c3-a3d4-4cf7-bee8-d1cfdac27177.png)


![15k_cf](https://user-images.githubusercontent.com/49332395/163659330-fa6c36ee-5367-4960-84a7-c0f3de4794d6.jpg)

# Usage Instructions
Make sure all the files in the src directory are in the same folder on your local machine. Open and run src/PR08D_V2.m, changing parameter as desired. The progam will automatically print all the output tables and show all of the graphs. 
