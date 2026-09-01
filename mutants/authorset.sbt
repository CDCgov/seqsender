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
            last "Lovelace",
            first "Ada",
            middle "Byron"
          }
        },
        {
          name name {
            last "Hopper",
            first "Grace",
            middle "Murray"
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
              last "Lovelace",
              first "Ada",
              middle "Byron"
            }
          },
          {
            name name {
              last "Hopper",
              first "Grace",
              middle "Murray"
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
      data str "Submission Title: sub1"
    }
  }
}
