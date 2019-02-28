%% List of all examples
example_dirs = {'1.3','1.5','2.5','2.6','6.1','6.2','6.3','6.4','8.1','9.1'};

%% Attempt to run all of the examples
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
fprintf('~~~ Running all of the examples ~~~\n')
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
for ii = 1:length(example_dirs)
    fprintf(['Test of exercise ',example_dirs{ii}, ' : '])
    try
        run(['../examples/exercise_',example_dirs{ii}, '/init.m']);
        fprintf('Pass\n')
    catch
        fprintf('Fail\n')
    end 
end