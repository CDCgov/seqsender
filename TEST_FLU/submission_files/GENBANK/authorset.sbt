Submit-block ::= {
  contact {
    contact {
      name name {
        last "Doe",
        first "Jane",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Centers for Disease Control and Prevention",
        div "Respiratory Viruses Branch, Division of Viral Diseases",
        city "Atlanta",
        sub "GA",
        country "USA",
        street "1600 Clifton Rd",
        email "email@myemail.com",
        phone "",
        postal-code "30329"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Doe",
            first "John"
          }
        },
        {
          name name {
            last "Doe",
            first "Jane"
          }
        }
      },
      affil std {
        affil "Centers for Disease Control and Prevention",
        div "Respiratory Viruses Branch, Division of Viral Diseases",
        city "Atlanta",
        sub "GA",
        country "USA",
        street "1600 Clifton Rd",
        postal-code "30329"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Doe",
              first "John"
            }
          },
          {
            name name {
              last "Doe",
              first "Jane"
            }
          }
        }
      },
      title "Test_Submission"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title: TEST_FLU"
    }
  }
}
