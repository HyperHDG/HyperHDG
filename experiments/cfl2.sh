python experiments/cfl2.py -i output/morgan-01-30.geo.h5 -b 100 --csv --cdf --noshow \
  | python experiments/plot.py -f csv -x ts -y cdf --log x \
               --xlabel 'timestep $t_s$' --tikz output/morgan-01-30.cfl.cdf
