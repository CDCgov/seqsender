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
                        'regex': '\s*[\S]+\s*'
                    },
                    'Password': {
                        'required': True,
                        'type': 'string',
                        'regex': '\s*[\S]+\s*'
                    },
                    'Spuid_Namespace': {
                        'required': True,
                        'type': 'string',
                        'regex': '\s*[\S]+\s*'
                    },
                    'BioSample_Package': {
                        'required': False,
                        'type': 'string',
                        'nullable': True
                    },
                    'GenBank_Auto_Remove_Failed_Samples': {
                        'required': False,
                        'type': 'boolean',
                        'nullable': True
                    },
                    'Publication_Title': {
                        'required': False,
                        'type': 'string',
                        'nullable': True
                    },
                    'Publication_Status': {
                        'required': False,
                        'type': 'string',
                        'regex': '(?i)(\W|^)(unpublished|in-press|published)(\W|$)',
                        'nullable': True
                    },
                    'Submission_Position': {
                        'required': False,
                        'type': 'integer',
                        'allowed': [1, 2],
                        'nullable': True
                    },
                    'Specified_Release_Date': {
                        'required': True,
                        'type': 'string',
                        'regex': '((?i)(\W|^)(\d+\s*(days|weeks|months)|\d{4}-\d{2}-\d{2})(\W|$))|(^\s*$)',
                        'nullable': True
                    },
                    'Link_Sample_Between_NCBI_Databases': {
                        'required': True,
                        'type': 'boolean',
                        'nullable': True
                    },
                    'Description': {
                        'required': True,
                        'type': 'dict',
                        'schema': {
                            'Organization': {
                                'required': True,
                                'type': 'dict',
                                'schema': {
                                    'Role': {
                                        'required': True,
                                        'type': 'string'
                                    },
                                    'Type': {
                                        'required': True,
                                        'type': 'string'
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
                                            'Postal_Code': {
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
                                                'required': False,
                                                'type': 'string',
                                                'nullable': True
                                            }
                                        }
                                    },
                                    'Submitter': {
                                        'required': True,
                                        'type': 'dict',
                                        'schema': {
                                            'Email': {
                                                'required': True,
                                                'type': 'string',
                                                'regex': '^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'
                                            },
                                            'Alt_Email': {
                                                'required': False,
                                                'type': 'string',
                                                'dependencies': ['Email'],
                                                'regex': '(^\s*$)|(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)',
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
                        'regex': '\s*[\S]+\s*',
                        'nullable': True
                    },
                    'Password': {
                        'required': False,
                        'type': 'string',
                        'regex': '\s*[\S]+\s*',
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
