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

% this function assigns a random value to the resistor
% given the tolerance rance. First generate a sign (either positive or negative),
% then generate a value in [r - r * tolerance, r + r * tolerance].
function [resistor] = make_random_resistor(r)
  s = rand;
  sign = NaN;
  if s < 0.5
    sign = -1;
  else
    sign = 1;
  end
  resistor = r + sign * (tolerance * s * r);
end

% using 0 based index for the resistor list, calculate v_0
function v_0 = calc_v_0(rs)
  n = rs(1) * rs(4);
  d = (rs(1) + rs(2)) * (rs(3) + rs(4)) + rs(1) * rs(2);
  v_0 = n / d;
end

function rs = assign_resistors()
    rs = [NaN,NaN,NaN,NaN];
    for i = 1:4
        rs(i) = make_random_resistor(resistors(i));
    end
end

% Calculate the voltage from the nominal resistor values
nominal_value = calc_v_0(resistors);

function plotHistogram(data)
    % Plot histogram of the data
    histogram(data);
    % Set title and axis labels
    title('Histogram of Data');
    xlabel('Data Values');
    ylabel('Frequency');
    % Save the histogram to a file in the current directory
    print('histogram100.png', '-dpng');
end


% Iterate through all the assignments and find the assignments producing the
% minimum and maximum voltage based on the randomly generated resistors
function [min, min_assignments, max, max_assignments, results] = compute_min_max(n)

    min = nominal_value;
    max = nominal_value;
    results = zeros(1,n);
    for i = 1:n
      rs = assign_resistors();
      v = calc_v_0(rs);
      results(i) = v;
      if v < min
        min = v;
        min_assignments = rs;
      elseif v > max
        max = v;
        max_assignments = rs;
      end
  end
end

[min, min_assignments, max, max_assignments, results] = compute_min_max(100);

nominal_value
min
min_assignments
max
max_assignments
plotHistogram(results);


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
