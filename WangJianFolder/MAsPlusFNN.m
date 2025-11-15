% Define the layers of the network
layers = [
    featureInputLayer(15009)  % Input layer with 15009 features per sample
    fullyConnectedLayer(10)    % Fully connected layer with 10 units
    reluLayer                  % ReLU activation function
    fullyConnectedLayer(26)    % Output layer with 26 units (assuming 26 classes)
    softmaxLayer               % Softmax activation function for multi-class classification
    classificationLayer];      % Classification layer for labels

% Specify Training Options
options = trainingOptions('sgdm', ...
    'MaxEpochs', 1000, ...
    'InitialLearnRate', 1e-3, ...
    'MiniBatchSize', 64, ...
    'Verbose', true, ...
    'Plots', 'training-progress');

% Example data loading (replace with your actual data loading method)
Data = load('Tumors_14.dat');  % Load your actual data file
XData = Data(:, 1:end-1);     % Features (input data)
YData = Data(:, end);         % Labels

% Normalize Data
XDataNormalized = (XData - min(XData)) ./ (max(XData) - min(XData));

% Example: Split data into training and test sets (adjust percentage as
% needed)
cv = cvpartition(size(XData, 1), 'HoldOut', 0.2); % 80% training, 20% test
idxTrain = cv.training;
idxTest = cv.test;

XTrain = XDataNormalized(idxTrain, :);
YTrain = YData(idxTrain);
XTest = XDataNormalized(idxTest, :);
YTest = YData(idxTest);

% Convert labels to categorical array for multi-class classification
YTrainCategorical = categorical(YTrain);

% Train the network
net = trainNetwork(XTrain, YTrainCategorical, layers, options);

% Evaluate the network on test data
YTestPred = classify(net, XTest);
YTest = categorical(YTest);
YTestPred = categorical(YTestPred);
accuracy = mean(YTestPred == YTest);
fprintf('Classification accuracy on test set: %.2f%%\n', accuracy * 100);