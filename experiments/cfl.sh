python experiments/cfl.py -i output/full.geo.h5 -b 100 --csv \
  | python experiments/plot.py -f csv -x ts -y density --log xy \
      --xlabel 'timestep $t_s$' --ylabel density --tikz output/cfl_hist
python experiments/cfl.py -i output/full.geo.h5 -b 100 --cdf --csv \
  | python experiments/plot.py -f csv -x ts -y cdf --log x \
      --xlabel 'timestep $t_s$' --ylabel 'cumulative density' --tikz output/cfl_cdf

