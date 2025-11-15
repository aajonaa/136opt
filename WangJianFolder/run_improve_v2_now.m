% Runner to execute ImproveTest_v2.m with a diary log, without modifying the script itself
try
    % Ensure we start from repo root
    repoRoot = fileparts(mfilename('fullpath')); % this file is in WangJianFolder
    repoRoot = fileparts(repoRoot); % go up to project root
    cd(repoRoot);

    % Start diary to a stable ASCII path relative to the repo
    logPath = fullfile('WangJianFolder', 'matlab_run_now.log');
    if exist(logPath, 'file')
        % Append separator for a new session
        dairy_fid = fopen(logPath, 'a'); %#ok<NASGU>
        fprintf('\n===== New session: %s =====\n', datestr(now));
        fclose('all');
    end
    diary(logPath); diary on;
    fprintf('RUN START: %s\n', datestr(now));

    % Execute the main script
    run('WangJianFolder/ImproveTest_v2.m');

    fprintf('RUN END: %s\n', datestr(now));
    diary off;
catch ME
    try diary on; catch, end
    disp(getReport(ME,'extended'));
    try diary off; catch, end
    rethrow(ME);
end

