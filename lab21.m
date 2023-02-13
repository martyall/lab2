% Problem 2.5

function main()


% define the resistors (in Ohms). Use 0 based indexing, so here
% r0 ~ R1 from the lab sheet.
r1 = 1000;
r2 = 2000;
r3 = 3000;
r4 = 4000;

resistors = [r1,r2,r3,r4];

% define the voltage source (in volts)
v_s = 10;

% each resistor has a 5% tolerance
tolerance = 0.05;

% using 0 based index for the resistor list, calculate v_0
function v_0 = calc_v_0(rs)
  n = rs(1) * rs(4);
  d = (rs(1) + rs(2)) * (rs(3) + rs(4)) + rs(1) * rs(2);
  v_0 = n / d;
end

% We take a subset of the resistors. If the resistor is in the subset, we set 
% it to its max resistence. If it's absent, we set it to its min.  
function rs = assign_resistors(max_set)
    rs = [NaN,NaN,NaN,NaN];
    for i = 1:4
      if ismember(i, max_set)
        rs(i) = resistors(i) + tolerance * resistors(i);
      else 
        rs(i) = resistors(i) - tolerance * resistors(i);
    end
end
end

% Calculate the voltage from the nominal resistor values
nominal_value = calc_v_0(resistors);

% Iterate through all the assignments and find the assignments producing the
% minimum and maximum voltage based on the +/-5% tolerance
function [min, min_assignments, max, max_assignments] = compute_min_max()

    min = nominal_value;
    max = nominal_value;
    for i = 1:4
      for combo = nchoosek([1,1,2,4],i)
        rs = assign_resistors(combo);
        v = calc_v_0(rs);
        if v < min
          min = v;
          min_assignments = rs;
        elseif v > max
          max = v;
          max_assignments = rs;
        end
    end
  end
end

[min, min_assignments, max, max_assignments] = compute_min_max();

nominal_value
min
min_assignments
max
max_assignments

% range
percent_range = [100 * ((min / nominal_value) - 1), 100 * ((max / nominal_value) - 1)]

end

%Solutions
%nominal_value =
%
%    0.1739
%
%
%min =
%
%    0.1620
%
%
%min_assignments =
%
%         950        2100        2850        3800
%
%
%max =
%
%    0.1846
%
%
%max_assignments =
%
%        1050        1900        2850        3800
%
%
%percent_range =
%
%   -6.8230    6.1538
