% Load dataset
Data = load('Tumors_14.dat');
XData = Data(:, 1:end-1);
YData = Data(:, end);

% Normalize Data
XDataNormalized = (XData - min(XData)) ./ (max(XData) - min(XData));

% Convert labels to categorical array for multi-class classification
YDataCategorical = categorical(YData);

% Ensure stratified split
cv = cvpartition(size(XData, 1), 'HoldOut', 0.2); % 80% training, 20% validation
idxTrain = cv.training;
idxValidation = cv.test;

% Ensure all classes are in the validation set
uniqueClasses = categories(YDataCategorical);
minSamplesPerClass = 2; % Adjust as needed

for i = 1:numel(uniqueClasses)
    class = uniqueClasses{i};
    classIndices = find(YDataCategorical == class);
    
    if sum(YDataCategorical(idxValidation) == class) < minSamplesPerClass
        % Move samples from training to validation
        numSamplesNeeded = minSamplesPerClass - sum(YDataCategorical(idxValidation) == class);
        validationCandidates = classIndices(idxTrain(classIndices));
        if numel(validationCandidates) < numSamplesNeeded
            error('Not enough samples to ensure every class is represented in the validation set.');
        end
        validationIndicesToAdd = validationCandidates(1:numSamplesNeeded);
        
        idxValidation(validationIndicesToAdd) = true;
        idxTrain(validationIndicesToAdd) = false;
    end
end

% Update training and validation sets
XTrain = XDataNormalized(idxTrain, :);
YTrainCategorical = YDataCategorical(idxTrain);
XValidation = XDataNormalized(idxValidation, :);
YValidation = YDataCategorical(idxValidation);

% Verify all classes are present in YValidation
classDistribution = tabulate(YValidation);

% Define the objective function to minimize validation loss
fobj = @(learningRate) trainAndEvaluateNetwork(XTrain, YTrainCategorical, XValidation, YValidation, learningRate);

% Run PSO to optimize learning rate
bestLearningRate = AO_Net(10, 50, 1, fobj);
% bestLearningRate = 1e-3;

% Train the network with the best learning rate
finalNet = trainNetwork(XTrain, YTrainCategorical, layers, ...
    trainingOptions('sgdm', ...
        'InitialLearnRate', bestLearningRate, ...
        'MaxEpochs', 100, ...
        'MiniBatchSize', 64, ...
        'Verbose', true, ...
        'Plots', 'none'));
    
% Evaluate final network on test set or further validate as needed
YTestPred = classify(finalNet, XValidation);
testAccuracy = mean(YTestPred == YValidation);
disp(['bestLearningRate:', num2str(bestLearningRate)]);
fprintf('Final network accuracy on validation set: %.2f%%\n', testAccuracy * 100);

% Evaluation function
function validationLoss = trainAndEvaluateNetwork(XTrain, YTrain, XValidation, YValidation, learningRate)
    % Define the layers of the network
    layers = [
        featureInputLayer(size(XTrain, 2))  % Input layer with correct number of features
        fullyConnectedLayer(10)    % Fully connected layer with 10 units
        reluLayer                  % ReLU activation function
        fullyConnectedLayer(26)    % Output layer with 26 units (assuming 26 classes)
        softmaxLayer               % Softmax activation function for multi-class classification
        classificationLayer        % Classification layer for labels
    ];

    % Specify Training Options
    options = trainingOptions('sgdm', ...
        'InitialLearnRate', learningRate, ...
        'MaxEpochs', 1000, ...
        'MiniBatchSize', 64, ...
        'Verbose', true, ...
        'Plots', 'none');

    % Train network
    net = trainNetwork(XTrain, YTrain, layers, options);

    % Validate the network
    YValidationPred = classify(net, XValidation);
    validationAccuracy = mean(YValidationPred == YValidation);

    % Convert categorical labels to one-hot encoding
    YValidationOneHot = onehotencode(YValidation, 2, 'single')'; % Ensure correct dimension

    % Get network predictions in probabilistic format
    YValidationPredProb = predict(net, XValidation)';

    % Ensure dimensions match
    if size(YValidationPredProb, 1) ~= size(YValidationOneHot, 1)
        error('The number of classes in predictions and actual labels do not match.');
    end

    % Calculate cross-entropy loss
    validationLoss = crossentropy(YValidationOneHot, YValidationPredProb);
    disp(['validation accuracy:', num2str(validationAccuracy * 100)]);
end