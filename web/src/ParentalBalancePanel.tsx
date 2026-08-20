/**
 * Whether chromosomal change falls on the two parental genomes differentially, or equally.
 *
 * The question the tool exists to answer, and the only panel that looks across the whole run
 * rather than at one sample. Everything else here reports events; this counts them.
 *
 * IT LEADS WITH WHAT IT CANNOT DO. Every path through `parentalBalance` that refuses returns a
 * headline saying so, and those are shown at the same size as a result, because the failure mode
 * this panel guards against is a reader taking "no difference found" for "the genomes are alike"
 * or taking a difference between unevenly-measured groups for biology.
 */
import { Text } from '@mantine/core'
import { parentalBalance, MIN_PER_GROUP, REPORTING_PER_GROUP } from './parentalBalance.ts'
import type { BalanceSample } from './parentalBalance.ts'

const pct = (x: number) => (Number.isFinite(x) ? `${(x * 100).toFixed(2)}%` : '-')

export function ParentalBalancePanel({ samples }: { samples: readonly BalanceSample[] }) {
  if (!samples.length) return null
  const r = parentalBalance(samples)
  const decided = r.verdict === 'differential' || r.verdict === 'equal'
  const colour = r.verdict === 'differential' ? 'var(--om-defect)' : 'var(--om-text-dim)'

  return (
    <div style={{
      border: '1px solid var(--om-line)', borderLeft: `4px solid ${colour}`,
      padding: '12px 15px', margin: '14px 0 4px',
    }}
    >
      <Text style={{ fontSize: 11, fontWeight: 700, letterSpacing: '0.06em', color: 'var(--om-text-dim)' }}>
        PARENTAL BALANCE, ACROSS THE WHOLE RUN
      </Text>
      <Text style={{ fontSize: 16, fontWeight: 700, lineHeight: 1.3, marginTop: 3, color: colour }}>
        {r.headline}
      </Text>

      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 18, marginTop: 10 }}>
        {r.groups.map((g) => (
          <div key={g.parent}>
            <Text style={{ fontSize: 11, fontWeight: 700, color: 'var(--om-text-dim)' }}>
              {g.parent.toUpperCase()} GENOMES
            </Text>
            <Text style={{ fontFamily: 'var(--om-mono)', fontSize: 13 }}>
              {g.samples} sample{g.samples === 1 ? '' : 's'}
            </Text>
            <Text style={{ fontFamily: 'var(--om-mono)', fontSize: 13 }}>
              {g.carryingAny} carrying any change
            </Text>
            <Text style={{ fontFamily: 'var(--om-mono)', fontSize: 13 }}>
              {g.events} distinct events, median {Number.isFinite(g.medianEvents) ? g.medianEvents : '-'}
            </Text>
            <Text size="xs" c="dimmed">
              noise ceiling {pct(g.medianNoiseCeiling)}
            </Text>
          </div>
        ))}
      </div>

      {/* The refusals are results. A reader who does not see these will read a null as a finding. */}
      <Text size="xs" c="dimmed" mt={9} style={{ maxWidth: 820, lineHeight: 1.55 }}>
        Groups matched on artefact propensity, not marker count: they differ by
        {' '}{Number.isFinite(r.power.skew) ? `${r.power.skew.toFixed(2)}x` : 'an unmeasurable amount'}
        {r.power.withinTolerance ? ', within tolerance' : ', past what this comparison allows'}.
        {' '}Exclusions {r.power.exclusion.balanced ? 'were balanced' : 'were NOT balanced'}
        {' '}({r.power.exclusion.maternal} maternal, {r.power.exclusion.paternal} paternal).
        {decided ? ` Smallest reachable p at these group sizes is ${r.minAchievableP.toFixed(4)}.` : ''}
        {' '}A hard floor of {MIN_PER_GROUP} genomes per group applies, and a result stays
        exploratory below {REPORTING_PER_GROUP}.
      </Text>

      {r.excluded.length > 0 && (
        <Text size="xs" c="dimmed" mt={6} style={{ maxWidth: 820, lineHeight: 1.55 }}>
          Excluded, and why: {r.excluded.map((e) => `${e.name} (${e.why})`).join('; ')}.
        </Text>
      )}

      {r.byClass.length > 0 && (
        <div style={{ marginTop: 10 }}>
          <Text style={{ fontSize: 11, fontWeight: 700, color: 'var(--om-text-dim)' }}>BY CLASS</Text>
          {r.byClass.map((c) => (
            <Text key={c.cls} style={{ fontFamily: 'var(--om-mono)', fontSize: 12, lineHeight: 1.7 }}>
              {c.cls}: {c.maternal} maternal against {c.paternal} paternal, p {c.p.toFixed(4)}
              {c.significant ? '  FLAGGED' : ''}
              {c.underpowered ? '  (could not have been flagged at this group size)' : ''}
            </Text>
          ))}
        </div>
      )}

      <Text size="xs" c="dimmed" mt={9} style={{ maxWidth: 820, lineHeight: 1.55 }}>
        {r.methods}
      </Text>
    </div>
  )
}
