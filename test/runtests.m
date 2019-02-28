%% List of all examples
example_dirs = {
    '1.3',...
    '1.5'};

%% Attempt to run all of the examples
for ii = 1:length(example_dirs)
    fprintf(['Test of exercise ',example_dirs{ii}, ' : '])
    try
        run(['../examples/exercise_',example_dirs{ii}, '/init.m']);
        fprintf('Pass\n')
    catch
        fprintf('Fail\n')
    end
    
end
