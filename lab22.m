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

% calculate v_0 using the specified values for the resistors
function v_0 = calc_v_0(rs)
  n = rs(1) * rs(4);
  d = (rs(1) + rs(2)) * (rs(3) + rs(4)) + rs(1) * rs(2);
  v_0 = n / d;
end

% this function assigns a random value to the resistor
% given the tolerance range. First generate a sign (either positive or negative),
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


% This function assigns a resistance value to each resistor 
% with the specified tolerance
function rs = assign_resistors()
    rs = [NaN,NaN,NaN,NaN];
    for i = 1:4
        rs(i) = make_random_resistor(resistors(i));
    end
end

% Calculate the voltage from the nominal resistor values
nominal_value = calc_v_0(resistors);

% Iterate through all the assignments and find the assignments producing the
% minimum and maximum voltage based on the randomly generated resistors
function [min, max, results, mean_value, std_deviation] = compute_min_max(n)

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

    % also return the mean and std_deviation of the distribution
    mean_value = mean(results);
    variance = var(results);
    std_deviation = sqrt(variance);
end
end

% This function plots a histogram of the data.
function plotHistogram(data)
    % Plot histogram of the data
    histogram(data);
    % Set title and axis labels
    htitle = sprintf('V0 histogram for %d trials', length(data));
    title(htitle);
    xlabel('V0 Values');
    ylabel('Frequency');
    % Save the histogram to a file in the current directory
    hname=sprintf('histogram%d.png', length(data));
    print(hname, '-dpng');
end

for m = [100,1000,10000,100000]

  [min, max, results, mean_value, std_deviation] = compute_min_max(m);
  m
  nominal_value
  min
  max
  mean_value
  std_deviation
  percent_range = [100 * ((min / nominal_value) - 1), 100 * ((max / nominal_value) - 1)]
  plotHistogram(results);
end
end
