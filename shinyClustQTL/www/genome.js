new Browser({
    chr:          'II',
    viewStart:    1000,
    viewEnd:      4000,
    cookieKey:    'yeast',

    coordSystem: {
      speciesName: 'S. cerevisiae',
      taxon: '9606',
      auth: '',
      version: '',
      ucscName: 'sacCer3'
    },

    sources:     [{name:                 'Genome',
                   twoBitURI:            'sacCer3.2bit',
                   tier_type:            'sequence',
                   provides_entrypoints: true},
                ]
  });
