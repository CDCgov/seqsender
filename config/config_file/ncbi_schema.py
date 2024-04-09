{
    'Submission': {
        'required': True,
        'type': 'dict',
        'schema': {
            'NCBI': {
                'required': True,
                'type': 'dict',
                'schema': {
                    'Username': {
                        'required': True,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+'
                    },
                    'Password': {
                        'required': True,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+'
                    },
                    'BioSample_Package': {
                        'required': False,
                        'type': 'string'
                    },
                    'Table2asn': {
                        'required': False,
                        'type': 'boolean'
                    },
                    'Submission_Position': {
                        'required': False,
                        'type': 'integer',
                        'allowed': [1, 2],
                        'nullable': True
                    },
                    'Description': {
                        'required': True,
                        'type': 'dict',
                        'schema': {
                            'Title': {
                                'required': True,
                                'type': 'string'
                            },
                            'Comment': {
                                'required': True,
                                'type': 'string'
                            },
                            'Organization': {
                                'required': True,
                                'type': 'dict',
                                'schema': {
                                    '@role': {
                                        'required': True,
                                        'type': 'string'
                                    },
                                    '@type': {
                                        'required': True,
                                        'type': 'string'
                                    },
                                    '@org_id': {
                                        'required': True,
                                        'type': 'integer'
                                    },
                                    'Name': {
                                        'required': True,
                                        'type': 'string'
                                    },
                                    'Address': {
                                        'required': True,
                                        'type': 'dict',
                                        'schema':{
                                            'Affil': {
                                                'required': True,
                                                'type': 'string'
                                            },
                                            'Div': {
                                                'required': True,
                                                'type': 'string'
                                            },
                                            'Street': {
                                                'required': True,
                                                'type': 'string'
                                            },
                                            'City': {
                                                'required': True,
                                                'type': 'string'
                                            },
                                            'Sub': {
                                                'required': True,
                                                'type': 'string'
                                            },
                                            'Postal_code': {
                                                'required': True,
                                                'type': 'integer'
                                            },
                                            'Country': {
                                                'required': True,
                                                'type': 'string'
                                            },
                                            'Email': {
                                                'required': True,
                                                'type': 'string',
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'
                                            },
                                            'Phone': {
                                                'required': True,
                                                'type': 'string',
                                                'nullable': True
                                            }
                                        }
                                    },
                                    'Submitter': {
                                        'required': True,
                                        'type': 'dict',
                                        'schema': {
                                            '@email': {
                                                'required': True,
                                                'type': 'string',
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'
                                            },
                                            '@alt_email': {
                                                'required': False,
                                                'type': 'string',
                                                'dependencies': ['@email'],
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$',
                                                'nullable': True
                                            },
                                            'Name': {
                                                'required': True,
                                                'type': 'dict',
                                                'schema': {
                                                    'First': {
                                                        'required': True,
                                                        'type': 'string'
                                                    },
                                                    'Last': {
                                                        'required': True,
                                                        'type': 'string'
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            },
            'GISAID': {
                'required': False,
                'type': 'dict',
                'schema': {
                    'Client-Id': {
                        'required': False,
                        'type': 'string',
                        'nullable': True
                    },
                    'Username': {
                        'required': False,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+',
                        'nullable': True
                    },
                    'Password': {
                        'required': False,
                        'type': 'string',
                        'regex': '^[a-zA-Z0-9_.+-]+',
                        'nullable': True
                    },
                    'Submission_Position': {
                        'required': False,
                        'type': 'integer',
                        'allowed': [1, 2],
                        'nullable': True
                    }
                }
            }
        }
    }
}
