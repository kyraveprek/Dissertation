function [start_row,end_row,info_row] = get_trial_rows(iT, pt, total_rows)
% get_trial_rows
%  gets the start row, end row, and info row from the subject data csvs
%  exported from cyberman or the quest

Trial1Start = 3;   %For the first trial for each subject the start row has to be set because the first two rows are 'info rows,' and have 'nans'
srowsubtract = 3;
erowsubtract = 1;
N_TRIALS = 122;

if iT == 1
    start_row = Trial1Start;
    end_row = pt(iT * 2) - 1;
    info_row = 1;
elseif iT == N_TRIALS
    start_row = pt(height(pt))+3;
    end_row = total_rows;
    info_row = pt((iT * 2) - 3) + 2;
else
    start_row = pt((iT * 2) - srowsubtract) + 4;
    end_row = pt((iT * 2) - erowsubtract) + 1;
    info_row = pt((iT * 2) - 3) + 2;
end

end