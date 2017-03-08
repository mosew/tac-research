function input = pb_7p2_example_u(s)
    input = 0.3*betapdf(s,12,7) + 0.6*betapdf(s,4,11);
end