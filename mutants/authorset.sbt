Submit-block ::= {
  contact {
    contact {
      name name {
        last "Lovelace",
        first "Ada"
      },
      affil std {
        affil "CDC",
        div "Division",
        city "Atlanta",
        sub "GA",
        country "USA",
        street "1 Main St",
        email "lab@example.org",
        postal-code "30329"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "",
            first "Lovelace"
          }
        },
        {
          name name {
            last "",
            first "Ada"
          }
        },
        {
          name name {
            last "",
            first "Hopper"
          }
        },
        {
          name name {
            last "",
            first "Grace"
          }
        }
      },
      affil std {
        affil "CDC",
        div "Division",
        city "Atlanta",
        sub "GA",
        country "USA",
        street "1 Main St",
        postal-code "30329"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "Unpublished",
      authors {
        names std {
          {
            name name {
              last "",
              first "Lovelace"
            }
          },
          {
            name name {
              last "",
              first "Ada"
            }
          },
          {
            name name {
              last "",
              first "Hopper"
            }
          },
          {
            name name {
              last "",
              first "Grace"
            }
          }
        }
      },
      title "Default publication title"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title: sub2"
    }
  }
}
