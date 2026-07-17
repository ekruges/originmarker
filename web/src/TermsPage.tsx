import { Anchor, Text, Title } from '@mantine/core'

const UPDATED = '2026-07-17'

export function TermsPage() {
  return (
    <div className="om-docs-wrap" style={{ maxWidth: 760, margin: '0 auto', padding: 12 }}>
      <article className="om-docs-body">
        <Title order={1} mb={4}>Terms of use</Title>
        <Text size="xs" c="dimmed" mb="lg">
          Last updated {UPDATED}. OriginMarker is free, open source, and offered as is.
        </Text>

        <h2 id="what">1. What this service is</h2>
        <p>
          OriginMarker proposes <strong>candidate</strong> flanking SNP markers around a
          pathogenic variant, assembled from public population data. Where primer design is
          enabled, it also proposes <strong>candidate</strong> PCR primer pairs for
          genotyping those markers. It is a research and educational tool.
        </p>
        <p>
          <strong>It is not a clinical diagnostic, and it is not medical advice.</strong> It
          does not diagnose, and it cannot determine which parental allele an embryo
          inherited. It cannot phase markers: that requires genotyping the family in a
          qualified genetics laboratory. Every marker it returns is a hypothesis to be
          validated, not a result.
        </p>
        <p>
          <strong>Primer pairs are candidates in the same sense.</strong> Nothing here runs a
          PCR. A pair is designed against a reference sequence, which is not any patient's
          genome: a variant private to a carrier can sit under a primer site, prevent that
          allele amplifying, and make a heterozygote read as a homozygote. The optional check
          against UCSC In-Silico PCR is an alignment against the same reference, not a
          wet-lab result, and a clean result from it is not a validation. Validate every pair
          at the bench before using it on a sample.
        </p>
        <p>
          Any clinical decision, including any decision about an embryo, is the
          responsibility of the qualified professionals making it. Nothing this tool outputs
          should be treated as a basis for such a decision on its own.
        </p>

        <h2 id="warranty">2. No warranty</h2>
        <p>
          The service and its source code are provided <strong>"as is", without warranty of
          any kind</strong>, express or implied, including but not limited to warranties of
          merchantability, fitness for a particular purpose, accuracy, and
          non-infringement. This mirrors the{' '}
          <Anchor href="https://www.apache.org/licenses/LICENSE-2.0" target="_blank" rel="noreferrer">
            Apache License 2.0
          </Anchor>{' '}
          the code is released under.
        </p>
        <p>
          The tool depends on third-party data that changes without notice and can be
          unavailable, incomplete, or wrong. Variant classifications are revised. Population
          frequencies are re-called. Genetic maps are population averages that describe no
          individual meiosis. Known limitations are documented in the{' '}
          <Anchor href="#/docs/not">Documentation</Anchor>, and that list is not exhaustive.
        </p>

        <h2 id="liability">3. Limitation of liability</h2>
        <p>
          To the maximum extent permitted by law, the author is not liable for any damages
          arising from use of this service or its output, including direct, indirect,
          incidental, special, exemplary, or consequential damages, and including any
          decision made in reliance on it.
        </p>

        <h2 id="use">4. Acceptable use</h2>
        <ul>
          <li>
            <strong>Do not submit patient-identifying information.</strong> The tool takes a
            variant identifier and nothing else. It has no field that needs a name, a date
            of birth, a medical record number, or any other identifier, and none should be
            entered anywhere, including the free-text box. This matters most in that box:
            text with no identifier in it is sent to a third-party model, verbatim, as
            section 5 describes.
          </li>
          <li>
            Do not present its output as a diagnostic result, or represent it as validated
            for clinical use. That includes primer pairs, which no part of this tool has run.
          </li>
          <li>
            Use it within the fair-use limits of the upstream data sources (below). Requests
            are rate limited. Do not attempt to circumvent those limits, and do not use the
            service in a way that degrades it for others.
          </li>
          <li>
            <strong>The primer check spends the operator's quota, not yours.</strong> UCSC
            publishes a limit of one request every 15 seconds and 5,000 a day for programmatic
            use, and this instance answers for that limit under its own key. Every pair you
            check, whether from the primer box or bundled into a build, spends part of it. It
            is rate limited per client for that reason. Run your own instance with your own
            key if you need it at any volume.
          </li>
        </ul>

        <h2 id="data">5. Third-party services, and what is sent to them</h2>
        <p>
          OriginMarker retrieves data at query time from{' '}
          <Anchor href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank" rel="noreferrer">ClinVar</Anchor>,{' '}
          <Anchor href="https://rest.ensembl.org" target="_blank" rel="noreferrer">Ensembl</Anchor>{' '}
          and{' '}
          <Anchor href="https://gnomad.broadinstitute.org" target="_blank" rel="noreferrer">gnomAD</Anchor>,
          and bundles the deCODE 2019 recombination map. Those datasets are governed by their
          own terms and licenses, which are not superseded by these terms or by the license
          on this code. Attribution and versions are recorded in every export and in the{' '}
          <Anchor href="#/docs/sources">Documentation</Anchor>.
        </p>
        <p>
          Two of those calls carry something you typed or something derived from it, so they
          are set out in full rather than left to the word "retrieves":
        </p>
        <ul>
          <li>
            <strong>Free text goes to a model, and only free text.</strong> If what you type
            contains an rsID, an HGVS expression or a ClinVar accession, a regular expression
            reads it here and <em>nothing is sent to any model</em>. If it contains no
            identifier, the text is sent verbatim to{' '}
            <Anchor href="https://www.anthropic.com" target="_blank" rel="noreferrer">Anthropic</Anchor>'s
            API, to be read by a small model that answers which variant you meant. That is a
            third party, on their terms and their retention policy, not this project's. It is
            why section 4 asks you to keep identifying information out of that box. Where the
            operator has configured no model key the box refuses free text instead, and the
            identifier paths keep working.
          </li>
          <li>
            <strong>Primer sequences go to UCSC, when you ask for the check.</strong> The
            optional verification sends the two primer sequences, the genome build, and the
            product size bounds to{' '}
            <Anchor href="https://genome.ucsc.edu/cgi-bin/hgPcr" target="_blank" rel="noreferrer">
              UCSC In-Silico PCR
            </Anchor>{' '}
            under this instance's API key. Those sequences are drawn from the reference
            genome, not from any sample, and no part of your query travels with them. Nothing
            is sent to UCSC unless you ask for the check.
          </li>
        </ul>
        <p>
          Those providers are not affiliated with OriginMarker and do not endorse it.
        </p>

        <h2 id="availability">6. Availability and privacy</h2>
        <p>
          This is a personal project running on personal hardware. There is no uptime
          commitment, no support commitment, and no guarantee that the service, its URLs, or
          its output format will continue to exist. Queries may be cached on the server to
          avoid re-requesting upstream data. No accounts, no tracking, no analytics.
        </p>
        <p>
          Your recent queries are kept in your own browser's local storage so the search box
          can offer them back. They are not sent anywhere, there is no account holding them,
          and Clear all removes them. Clearing your browser's site data removes them too.
        </p>
        <p>
          "No tracking" is about who you are, and is not a claim that nothing leaves the
          server. Section 5 lists what does: every build queries ClinVar, Ensembl and gnomAD
          for the variant you named; free text with no identifier in it is sent to a
          third-party model; primer sequences are sent to UCSC if you ask for the check. None
          of those carries anything about a patient unless you type it into the box, which is
          why section 4 asks you not to.
        </p>

        <h2 id="license">7. License</h2>
        <p>
          The source code is released under the{' '}
          <Anchor href="https://github.com/ekruges/originmarker/blob/main/LICENSE" target="_blank" rel="noreferrer">
            Apache License 2.0
          </Anchor>. You may use, modify and redistribute it, including commercially, subject
          to that license. Running your own instance is encouraged and is the intended way to
          use it at any scale.
        </p>
        <p>
          <strong>One dependency is not Apache 2.0, and it is optional for that reason.</strong>{' '}
          Primer design needs{' '}
          <Anchor href="https://github.com/libnano/primer3-py" target="_blank" rel="noreferrer">primer3-py</Anchor>,
          which is GPLv2. It is deliberately absent from <code>requirements.txt</code> and from
          the default image, so what this repository distributes is Apache 2.0 throughout.
          Installing it on a server you run is not distribution and triggers nothing. If you
          build an image with it switched on and then redistribute that image, the combined
          work is subject to the GPLv2, and that is your decision to make rather than this
          project's default. Nothing from UCSC's own isPcr source is included here: the
          optional check calls their hosted service.
        </p>

        <h2 id="changes">8. Changes</h2>
        <p>
          These terms may change. The date at the top reflects the last revision, and the
          history is public in the{' '}
          <Anchor href="https://github.com/ekruges/originmarker" target="_blank" rel="noreferrer">
            repository
          </Anchor>.
        </p>

        <h2 id="contact">9. Contact</h2>
        <p>
          <Anchor href="mailto:kruger.ezra.s@gmail.com">kruger.ezra.s@gmail.com</Anchor>
        </p>
      </article>
    </div>
  )
}
