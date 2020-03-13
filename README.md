# Hand Flipping classification code 

# Instructions to run the flipping analysis code


 1. Install latest version of Matlab
 2. Clone this github repository on your machine or download project as ZIP
 3. Open the matlab code ["hand_flip_bookchap.m"](https://github.com/wearablebiosensing/wearablesensorsbook/blob/master/Matlab%20Code/hand_flip_bookchap.m "hand_flip_bookchap.m") in the folder Matlab Code
 4. Set the variable "analysis_side" as 0 for Right hand analysis or 1 for Left hand
 5. Run the code and observe that a new file titled "features_handflip_*_labelled_10subs.mat" is saved where * refers to analysis side
 
# Instructions to perform classification
1. Load the "features_handflip_*_labelled_10subs.mat" using [Matlab Load](https://www.mathworks.com/help/matlab/ref/load.html) command
2. Open Classification Learner from the Apps section of Matlab
![image](https://user-images.githubusercontent.com/1295467/76643949-72f75680-652c-11ea-8249-ed34537b33e9.png)

3. Create New Session with workspace
![2020-03-13](https://user-images.githubusercontent.com/1295467/76644258-03359b80-652d-11ea-8891-188eb2fed661.png)

4. Import features variables into classification learner, Set response to Label and ensure the subject_no checkbox is unselected in the list of predictors
![image](https://user-images.githubusercontent.com/1295467/76644446-5b6c9d80-652d-11ea-8ccd-a60eeafe929c.png)

5. Create a session and select Fine Gaussian SVM 
![image](https://user-images.githubusercontent.com/1295467/76644834-2ad93380-652e-11ea-9392-d0ccde32e4b8.png)
6. Click on train button to train the classifier
