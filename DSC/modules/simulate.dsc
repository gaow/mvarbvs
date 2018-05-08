
original_Y: Python(data['Y'] = numpy.vstack(data['Y'].values()).T)
  # do not simulate data, just use original
  data: $data
  $data: data
