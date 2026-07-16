import { Anchor, Text, Title } from '@mantine/core'

const UPDATED = '2026-07-16'

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
          pathogenic variant, assembled from public population data. It is a research and
          educational tool.
        </p>
        <p>
          <strong>It is not a clinical diagnostic, and it is not medical advice.</strong> It
          does not diagnose, and it cannot determine which parental allele an embryo
          inherited. It cannot phase markers: that requires genotyping the family in a
          qualified genetics laboratory. Every marker it returns is a hypothesis to be
          validated, not a result.
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
            entered anywhere, including the free-text box.
          </li>
          <li>
            Do not present its output as a diagnostic result, or represent it as validated
            for clinical use.
          </li>
          <li>
            Use it within the fair-use limits of the upstream data sources (below). Requests
            are rate limited. Do not attempt to circumvent those limits, and do not use the
            service in a way that degrades it for others.
          </li>
        </ul>

        <h2 id="data">5. Third-party data</h2>
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
          Those providers are not affiliated with OriginMarker and do not endorse it.
        </p>

        <h2 id="availability">6. Availability and privacy</h2>
        <p>
          This is a personal project running on personal hardware. There is no uptime
          commitment, no support commitment, and no guarantee that the service, its URLs, or
          its output format will continue to exist. Queries may be cached on the server to
          avoid re-requesting upstream data. No accounts, no tracking, no analytics.
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
