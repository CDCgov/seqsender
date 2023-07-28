Submit-block ::= {
  contact {
    contact {
      name name {
        last "Doe",
        first "John",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "NIH",
        div "NCBI",
        city "Bethesda",
        sub "MD",
        country "USA",
        street "10 Center Dr",
        email "jdoe@nih.gov",
        phone "301-402-8219",
        postal-code "20895"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Doe",
            first "John",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "NIH",
        div "NCBI",
        city "Bethesda",
        sub "MD",
        country "USA",
        street "10 Center Dr",
        postal-code "20895"
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
              first "John",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "SARS-CoV2 sequence detection"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL: jdoe@nih.gov"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title: None"
    }
  }
}
